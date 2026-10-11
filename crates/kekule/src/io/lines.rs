//! Line reading shared by streaming record readers.

use std::io::BufRead;

/// Reads one line ended by `\n`, `\r\n`, or `\r` into `line`, keeping at most
/// `limit` bytes. Returns `None` at the end of input, or whether bytes were
/// dropped.
pub(crate) fn read_line(
    reader: &mut impl BufRead,
    line: &mut Vec<u8>,
    limit: usize,
) -> std::io::Result<Option<bool>> {
    line.clear();
    let mut truncated = false;
    let mut read_any = false;
    loop {
        let (consumed, terminator) = {
            let buffer = match reader.fill_buf() {
                Ok(buffer) => buffer,
                Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(error) => return Err(error),
            };
            if buffer.is_empty() {
                return Ok(read_any.then_some(truncated));
            }
            read_any = true;
            let end = buffer
                .iter()
                .position(|&byte| byte == b'\n' || byte == b'\r');
            let content = &buffer[..end.unwrap_or(buffer.len())];
            let room = limit.saturating_sub(line.len());
            truncated |= content.len() > room;
            line.extend_from_slice(&content[..content.len().min(room)]);
            match end {
                Some(end) => (end + 1, Some(buffer[end])),
                None => (buffer.len(), None),
            }
        };
        reader.consume(consumed);
        match terminator {
            None => continue,
            Some(b'\r') => loop {
                // `\r\n` is one terminator.
                match reader.fill_buf() {
                    Ok(next) => {
                        if next.first() == Some(&b'\n') {
                            reader.consume(1);
                        }
                        break;
                    }
                    Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
                    Err(error) => return Err(error),
                }
            },
            Some(_) => {}
        }
        return Ok(Some(truncated));
    }
}
