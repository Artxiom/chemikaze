use crate::atom_counts::AtomCounts;
use crate::errors::ChemikazeError;
use crate::errors::ErrorKind::Parsing;
use crate::periodic_table;
use crate::periodic_table::EARTH_ELEMENT_CNT;

const OP: u8 = '(' as u8;
const CP: u8 = ')' as u8;
const DOT: u8 = '.' as u8;
const ZERO: u8 = '0' as u8;

const MF_PUNCTUATION: [u8; 7] = [OP, CP, '+' as u8, '-' as u8, DOT, '[' as u8, ']' as u8];

const DEFAULT_BUFFER_SIZE: usize = 128;

/// ParseInput represents a molecular formula input.
pub struct ParseInput<'a>(&'a [u8]);

/// SanitizedInput represents a sanitized molecular formula input (surrounding whitespaces removed, non-empty)
pub struct SanitizedInput<'a>(&'a [u8]);

impl<'a> From<&'a [u8]> for ParseInput<'a> {
    fn from(input: &'a [u8]) -> Self {
        Self(input)
    }
}

impl<'a> From<&'a str> for ParseInput<'a> {
    fn from(input: &'a str) -> Self {
        Self(input.as_bytes())
    }
}

/// Creates a new `SanitizedInput` from a `ParseInput`.
impl<'a> TryFrom<ParseInput<'a>> for SanitizedInput<'a> {
    type Error = ChemikazeError;

    fn try_from(input: ParseInput<'a>) -> Result<Self, Self::Error> {
        let input = input.0.trim_ascii();
        if input.is_empty() {
            return Err(ChemikazeError {
                kind: Parsing,
                msg: "Empty Molecular Formula".into(),
            });
        }

        Ok(Self(input))
    }
}

/// Molecular Formula Parser
pub struct MfParser {
    pos: usize,
    elements: Vec<u8>,
    coeffs: Vec<u32>,
}

impl MfParser {
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            pos: 0,
            elements: Vec::with_capacity(capacity),
            coeffs: Vec::with_capacity(capacity),
        }
    }

    pub fn new() -> Self {
        Self::with_capacity(DEFAULT_BUFFER_SIZE)
    }

    #[allow(dead_code)]
    pub fn parse_single<'a>(
        input: impl Into<ParseInput<'a>>,
    ) -> Result<AtomCounts, ChemikazeError> {
        let input = input.into();
        let sanitized_input = SanitizedInput::try_from(input)?;
        let mut parser = MfParser::with_capacity(sanitized_input.0.len());
        parser.parse_sanitized(sanitized_input)
    }

    pub fn parse<'a>(
        &mut self,
        input: impl Into<ParseInput<'a>>,
    ) -> Result<AtomCounts, ChemikazeError> {
        let input = SanitizedInput::try_from(input.into())?;
        self.parse_sanitized(input)
    }

    fn parse_sanitized(&mut self, input: SanitizedInput) -> Result<AtomCounts, ChemikazeError> {
        let input = input.0;

        // Initialize working vectors
        self.coeffs.clear();
        self.elements.clear();
        self.elements.resize(input.len(), 0);
        self.coeffs.resize(input.len(), 0);

        // First pass: read symbols and their immediate coefficients
        self.pos = 0;
        self.read_symbols_and_coeffs(input)?;

        // Second pass: apply group coefficients (parentheses and leading numbers)
        self.pos = 0;
        self.apply_group_coefficients(input)?;

        // Combine into final counts
        let mut counts = [0u32; EARTH_ELEMENT_CNT];
        self.elements
            .iter()
            .zip(self.coeffs.iter())
            .for_each(|(&element, &coeff)| {
                if coeff > 0 {
                    counts[element as usize] += coeff;
                }
            });

        Ok(AtomCounts { counts })
    }

    fn read_symbols_and_coeffs(&mut self, input: &[u8]) -> Result<(), ChemikazeError> {
        while self.pos < input.len() {
            let letter = input[self.pos];

            if letter.is_ascii_uppercase() {
                let start_pos = self.pos;
                let element = self.parse_element_symbol(input)?;
                let coeff = self.parse_coefficient(input);

                self.elements[start_pos] = element;
                self.coeffs[start_pos] = coeff;
            } else if letter.is_ascii_digit() || MF_PUNCTUATION.contains(&letter) {
                self.pos += 1;
            } else {
                return Err(ChemikazeError {
                    kind: Parsing,
                    msg: format!("Unexpected symbol: '{}'", char::from(letter)),
                });
            }
        }
        Ok(())
    }

    fn parse_element_symbol(&mut self, input: &[u8]) -> Result<u8, ChemikazeError> {
        if self.pos >= input.len() {
            return Err(ChemikazeError {
                kind: Parsing,
                msg: "Unexpected end of input".into(),
            });
        }

        let first = input[self.pos];
        self.pos += 1;

        let symbol = if self.pos < input.len() && input[self.pos].is_ascii_lowercase() {
            let second = input[self.pos];
            self.pos += 1;
            [first, second]
        } else {
            [first, 0]
        };

        periodic_table::get_element_by_symbol_bytes(symbol)
    }

    fn parse_coefficient(&mut self, input: &[u8]) -> u32 {
        if self.pos >= input.len() || !input[self.pos].is_ascii_digit() {
            return 1;
        }

        let mut coeff = 0u32;
        while self.pos < input.len() && input[self.pos].is_ascii_digit() {
            coeff = coeff * 10 + (input[self.pos] - ZERO) as u32;
            self.pos += 1;
        }
        coeff
    }

    fn apply_group_coefficients(&mut self, input: &[u8]) -> Result<(), ChemikazeError> {
        self.pos = 0;
        let mut paren_depth = 0i32;

        while self.pos < input.len() {
            // Handle leading coefficients like "2H2O"
            if input[self.pos].is_ascii_digit() {
                let coeff = self.parse_coefficient(input);
                self.scale_forward(input, coeff, paren_depth);
            }

            // Skip alphanumeric characters
            while self.pos < input.len() && input[self.pos].is_ascii_alphanumeric() {
                self.pos += 1;
            }

            if self.pos >= input.len() {
                break;
            }

            match input[self.pos] {
                OP => {
                    paren_depth += 1;
                    self.pos += 1;
                }
                CP => {
                    let close_pos = self.pos;
                    self.pos += 1;
                    let coeff = self.parse_coefficient(input);
                    self.scale_backward(input, close_pos, coeff, paren_depth);
                    paren_depth -= 1;
                }
                _ => self.pos += 1,
            }
        }

        if paren_depth != 0 {
            return Err(ChemikazeError {
                kind: Parsing,
                msg: "The opening and closing parentheses don't match.".into(),
            });
        }

        Ok(())
    }

    fn scale_forward(&mut self, input: &[u8], group_coeff: u32, curr_depth: i32) {
        if group_coeff == 1 {
            return;
        }

        let mut depth = curr_depth;
        let mut pos = self.pos;

        while pos < input.len() && depth >= curr_depth {
            match input[pos] {
                OP => depth += 1,
                CP => depth -= 1,
                DOT if depth == curr_depth => break,
                _ => {}
            }

            self.coeffs[pos] = self.coeffs[pos].saturating_mul(group_coeff);
            pos += 1;
        }
    }

    fn scale_backward(
        &mut self,
        input: &[u8],
        close_pos: usize,
        group_coeff: u32,
        curr_depth: i32,
    ) {
        if group_coeff == 1 {
            return;
        }

        let mut depth = curr_depth;
        let mut pos = close_pos;

        while pos > 0 && depth <= curr_depth {
            pos -= 1;

            match input[pos] {
                OP => depth += 1,
                CP => depth -= 1,
                _ => {}
            }

            self.coeffs[pos] = self.coeffs[pos].saturating_mul(group_coeff);
        }
    }
}

#[cfg(test)]
mod parse_mf_test {
    use super::*;

    //#[inline(always)]
    pub fn parse_mf<'a>(input: impl Into<ParseInput<'a>>) -> Result<AtomCounts, ChemikazeError> {
        MfParser::parse_single(input)
    }

    #[test]
    fn simple_mf_is_parsed_into_counts() {
        assert_eq!("H2O", parse_mf("H2O").unwrap().to_string());
        assert_eq!("H2O", parse_mf("HOH").unwrap().to_string());
        assert_eq!("H132C67O3N8", parse_mf("C67H132N8O3").unwrap().to_string());
    }

    #[test]
    fn complicated_mf_is_parsed_into_counts() {
        assert_eq!(
            "H12O6NSCl3Na3",
            parse_mf("[(2H2O.NaCl)3S.N]2-").unwrap().to_string()
        );
        assert_eq!(
            "H12O6NSCl3Na3",
            parse_mf(" [(2H2O.NaCl)3S.N]2- ").unwrap().to_string()
        );
    }

    #[test]
    fn errs_on_empty_mf() {
        assert_eq!("Empty Molecular Formula", parse_mf("").unwrap_err().msg);
        assert_eq!("Empty Molecular Formula", parse_mf(" ").unwrap_err().msg);
        assert_eq!("Empty Molecular Formula", parse_mf("  ").unwrap_err().msg);
    }

    #[test]
    fn trims_input() {
        assert_eq!("H8C2", parse_mf("  CH4CH4 ").unwrap().to_string());
        assert_eq!("H5C2", parse_mf("  (CH4).[CH]-  ").unwrap().to_string());
    }

    #[test]
    fn parenthesis_multiply_counts() {
        assert_eq!("H8C2", parse_mf("(CH4CH4)").unwrap().to_string());
        assert_eq!("H16C4", parse_mf("(CH4CH4)2").unwrap().to_string());
        assert_eq!("H16C5", parse_mf("C(CH4CH4)2").unwrap().to_string());
        assert_eq!("H4C2O4P", parse_mf("(C(OH)2)2P").unwrap().to_string());
        assert_eq!("C2O2PS8", parse_mf("(C(2S)2O)2P").unwrap().to_string());
        assert_eq!(
            "H2C2O2PS4",
            parse_mf("(C(OH))2(S(S))2P").unwrap().to_string()
        );
    }

    #[test]
    fn number_at_the_beginning_multiples_counts() {
        assert_eq!("H4O2", parse_mf("2H2O").unwrap().to_string());
        assert_eq!("", parse_mf("0H2O").unwrap().to_string());
    }

    #[test]
    fn sign_is_ignored_in_counts() {
        assert_eq!("H8C2", parse_mf("[CH4CH4]+").unwrap().to_string());
        assert_eq!("H8C2", parse_mf("[CH4CH4]2+").unwrap().to_string());
    }

    #[test]
    fn dots_separate_components_but_components_are_summed_up() {
        assert_eq!("H6CN", parse_mf("NH3.CH3").unwrap().to_string());
        assert_eq!("H9C2N", parse_mf("NH3.2CH3").unwrap().to_string());
        assert_eq!("H9C2N", parse_mf("2CH3.NH3").unwrap().to_string());
    }

    #[test]
    fn errs_if_parenthesis_do_not_match() {
        assert_eq!(
            "The opening and closing parentheses don't match.",
            parse_mf("(C").unwrap_err().msg
        );
        assert_eq!(
            "The opening and closing parentheses don't match.",
            parse_mf(")C").unwrap_err().msg
        );
        assert_eq!(
            "The opening and closing parentheses don't match.",
            parse_mf("C)").unwrap_err().msg
        );
        assert_eq!(
            "The opening and closing parentheses don't match.",
            parse_mf("C(").unwrap_err().msg
        );
        assert_eq!(
            "The opening and closing parentheses don't match.",
            parse_mf("(C))").unwrap_err().msg
        );
        assert_eq!(
            "The opening and closing parentheses don't match.",
            parse_mf("(C(OH)2(S(S))2P").unwrap_err().msg
        );
    }

    #[test]
    fn errs_if_symbol_is_invalid() {
        use crate::errors::ErrorKind::UnknownElement;
        let err = parse_mf("A").err().unwrap();
        assert_eq!(UnknownElement, err.kind);
        assert!(err.msg.contains("Unknown chemical symbol: A"));

        let err = parse_mf("o").err().unwrap();
        assert_eq!(Parsing, err.kind);
        assert_eq!("Unexpected symbol: 'o'", err.msg);
    }

    #[test]
    fn errs_if_contains_special_symbols_outside_of_allowed_punctuation() {
        let err = parse_mf("=").err().unwrap();
        assert_eq!(Parsing, err.kind);

        let err = parse_mf("O=").unwrap_err();
        assert_eq!(Parsing, err.kind);

        let err = parse_mf("=C").unwrap_err();
        assert_eq!(Parsing, err.kind);
    }
}
