use std::{cmp, fmt};
use std::iter::FromIterator;

#[derive(Debug, Copy, Clone, Hash, Eq, PartialEq)]
pub struct Range {
    pub start: u64,
    pub end: u64,
}

impl Range {
    pub fn new(start: u64, end: u64) -> Range {
        Range {
            start: start,
            end: end,
        }
    }

    pub fn overlaps(&self, other: &Range) -> bool {
        (other.start >= self.start && other.start <= self.end)
        || (other.end >= self.start && other.end <= self.end)
    }

    pub fn contains(&self, point: u64) -> bool {
        point >= self.start && point <= self.end
    }

    pub fn merge(&mut self, other: &Range) {
        self.start = cmp::min(self.start, other.start);
        self.end = cmp::max(self.end, other.end);
    }
}

impl fmt::Display for Range {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[{},{}]", self.start, self.end)
    }
}

#[derive(Debug, Clone)]
pub struct RangeStack {
    pub ranges: Vec<Range>,
}

impl RangeStack {
    fn add(&mut self, range: &Range) {
        if let Some(last) = self.ranges.last_mut() {
            if last.overlaps(range) {
                last.merge(range);
                return;
            }
        }

        self.ranges.push(*range);
    }
}


impl FromIterator<Range> for RangeStack {
    fn from_iter<I>(iterator: I) -> Self
        where I: IntoIterator<Item=Range>
    {
        let mut raw_ranges: Vec<_> = iterator.into_iter().collect();
        raw_ranges.sort_by(|a,b| a.start.cmp(&b.start));

        let mut range_stack = RangeStack {
            ranges: Vec::new(),
        };

        for range in &raw_ranges {
            range_stack.add(range);
        }

        range_stack
    }
}

impl<'a> FromIterator<&'a Range> for RangeStack {
    fn from_iter<I>(iterator: I) -> Self
        where I: IntoIterator<Item=&'a Range>
    {
        iterator.into_iter().cloned().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::Range;

    #[test]
    fn test_range_contains() {
        let range = Range::new(100,110);
        assert_eq!(range.contains(100), true);
        assert_eq!(range.contains(105), true);
        assert_eq!(range.contains(110), true);
        assert_eq!(range.contains(99), false);
        assert_eq!(range.contains(111), false);
    }

    #[test]
    fn test_range_display() {
        let range = Range::new(100,110);
        assert_eq!(format!("{}", range), "[100,110]");
    }
}