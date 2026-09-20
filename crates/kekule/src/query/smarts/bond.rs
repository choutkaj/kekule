use super::*;

pub(super) fn parse(
    input: &str,
    start: usize,
    options: SmartsParseOptions,
) -> Result<(BondExpression, Option<bool>, usize), SmartsParseError> {
    let mut parser = BondParser {
        input: input.as_bytes(),
        cursor: start,
        direction: None,
    };
    let expression = parser.low_and()?;
    if expression.node_count() > options.max_expression_nodes
        || expression.depth() > options.max_expression_depth
    {
        return Err(SmartsParseError::limit(
            start..parser.cursor,
            "bond expression nodes",
            expression.node_count(),
            options.max_expression_nodes,
        ));
    }
    Ok((expression, parser.direction, parser.cursor))
}

struct BondParser<'a> {
    input: &'a [u8],
    cursor: usize,
    direction: Option<bool>,
}
impl BondParser<'_> {
    fn peek(&self) -> Option<u8> {
        self.input.get(self.cursor).copied()
    }
    fn low_and(&mut self) -> Result<BondExpression, SmartsParseError> {
        let mut terms = vec![self.or()?];
        while self.peek() == Some(b';') {
            self.cursor += 1;
            terms.push(self.or()?);
        }
        BondExpression::all(terms).map_err(|e| expression_error(0..self.cursor, e))
    }
    fn or(&mut self) -> Result<BondExpression, SmartsParseError> {
        let mut terms = vec![self.and()?];
        while self.peek() == Some(b',') {
            self.cursor += 1;
            terms.push(self.and()?);
        }
        BondExpression::any(terms).map_err(|e| expression_error(0..self.cursor, e))
    }
    fn and(&mut self) -> Result<BondExpression, SmartsParseError> {
        let mut terms = vec![self.unary()?];
        loop {
            if self.peek() == Some(b'&') {
                self.cursor += 1;
            } else if !self.peek().is_some_and(is_bond_start) {
                break;
            }
            terms.push(self.unary()?);
        }
        BondExpression::all(terms).map_err(|e| expression_error(0..self.cursor, e))
    }
    fn unary(&mut self) -> Result<BondExpression, SmartsParseError> {
        let start = self.cursor;
        let mut negate = false;
        while self.peek() == Some(b'!') {
            self.cursor += 1;
            negate = !negate;
        }
        let byte = self.peek().ok_or_else(|| {
            SmartsParseError::syntax(start..self.cursor, "missing bond primitive")
        })?;
        self.cursor += 1;
        let expression = match byte {
            b'-' => non_aromatic_bond_expression(BondOrder::Single),
            b'=' => non_aromatic_bond_expression(BondOrder::Double),
            b'#' => Ok(BondExpression::predicate(BondPredicate::Order(
                BondOrder::Triple,
            ))),
            b'$' => Ok(BondExpression::predicate(BondPredicate::Order(
                BondOrder::Quadruple,
            ))),
            b':' => aromatic_bond_expression(),
            b'~' => Ok(BondExpression::always()),
            b'@' => Ok(BondExpression::predicate(BondPredicate::RingMembership(
                true,
            ))),
            b'/' | b'\\' => {
                // Directional alternatives are normalized by the stereo parser.
                let up = byte == b'/';
                if self.peek() == Some(b'?') {
                    return Err(SmartsParseError::unsupported(
                        start..self.cursor + 1,
                        "unspecified directional syntax is outside this dialect",
                    ));
                }
                if self.direction.is_none() {
                    self.direction = Some(up);
                }
                BondExpression::all([
                    default_bond_expression()?,
                    BondExpression::predicate(BondPredicate::Direction(up)),
                ])
            }
            _ => {
                return Err(SmartsParseError::syntax(
                    start..self.cursor,
                    "missing or invalid bond primitive",
                ))
            }
        }
        .map_err(|e| expression_error(start..self.cursor, e))?;
        if negate {
            expression
                .negate()
                .map_err(|e| expression_error(start..self.cursor, e))
        } else {
            Ok(expression)
        }
    }
}
