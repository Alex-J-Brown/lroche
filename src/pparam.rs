use std::str::FromStr;
use serde::{Serialize, Deserialize};


#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Pparam {
    pub value: f64,
    pub range: f64,
    pub dstep: f64,
    pub vary: bool,
    pub defined: bool
}


impl Default for Pparam {
    fn default() -> Self {
        Self {
            value: 0.0,
            range: 0.0,
            dstep: 0.0,
            vary: false,
            defined: false
        }
    }
}


impl FromStr for Pparam {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split_whitespace();
        let _name: &str = fields.next().ok_or("missing value")?;
        let _equals: &str = fields.next().ok_or("missing value")?;
        let value = fields.next().ok_or("missing value")?.parse().map_err(|_| "bad value")?;
        let range = fields.next().ok_or("missing range")?.parse().map_err(|_| "bad range")?;
        let dstep = fields.next().ok_or("missing dstep")?.parse().map_err(|_| "bad dstep")?;
        let vary = fields.next().ok_or("missing vary")?.parse::<i32>().map_err(|_| "bad vary")? != 0;
        let defined = fields.next().ok_or("missing vary")?.parse::<i32>().map_err(|_| "bad vary")? != 0;

        Ok(Pparam { value, range, dstep, vary: vary, defined: defined})
    }
}
