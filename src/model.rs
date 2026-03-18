use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::{f64::consts::PI, str::FromStr};
use crate::vec3::Vec3;


#[derive(Debug)]
pub enum Entry {
    Param(Pparam),
    Scalar(String),
}


fn parse_entry(line: &str) -> Option<(String, Entry)> {
    let line = line.trim();
    if line.is_empty() {
        return None;
    }

    let mut it = line.split_whitespace();

    let name = it.next()?.to_string();
    if it.next()? != "=" {
        return None;
    }

    let rest: Vec<_> = it.collect();

    match rest.len() {
        1 => Some((name, Entry::Scalar(rest[0].to_string()))),
        5 => {
            let joined = format!("{} = {}", name, rest.join(" "));
            let p: Pparam = joined.parse().ok()?;
            Some((name, Entry::Param(p)))
        }
        _ => None,
    }
}

fn read_lines<P: AsRef<Path>>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>> {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn load_entries(path: &str) -> Result<HashMap<String, Entry>, String> {
    let mut map = HashMap::new();

    for line in read_lines(path).map_err(|e| e.to_string())?.map_while(Result::ok) {
        if let Some((name, entry)) = parse_entry(&line) {
            map.insert(name, entry);
        }
    }

    Ok(map)
}

fn get_p(map: &HashMap<String, Entry>, k: &str) -> Result<Pparam, String> {
    match map.get(k) {
        Some(Entry::Param(p)) => Ok(*p),
        _ => Err(format!("missing Pparam: {}", k)),
    }
}

fn get_f64(map: &HashMap<String, Entry>, k: &str) -> Result<f64, String> {
    match map.get(k) {
        Some(Entry::Scalar(v)) => v.parse().map_err(|_| format!("bad f64: {}", k)),
        _ => Err(format!("missing f64: {}", k)),
    }
}

fn get_u32(map: &HashMap<String, Entry>, k: &str) -> Result<u32, String> {
    match map.get(k) {
        Some(Entry::Scalar(v)) => v.parse().map_err(|_| format!("bad u32: {}", k)),
        _ => Err(format!("missing u32: {}", k)),
    }
}

fn get_bool(map: &HashMap<String, Entry>, k: &str) -> Result<bool, String> {
    match map.get(k) {
        Some(Entry::Scalar(v)) => match v.as_str() {
            "0" => Ok(false),
            "1" => Ok(true),
            _ => Err(format!("bad bool: {}", k)),
        },
        _ => Err(format!("missing bool: {}", k)),
    }
}

fn get_ldc(map: &HashMap<String, Entry>, k: &str) -> Result<LDCType, String> {
    match map.get(k) {
        Some(Entry::Scalar(v)) => match v.as_str() {
            "Claret" => Ok(LDCType::Claret),
            "Poly" => Ok(LDCType::Poly),
            _ => Err(format!("bad LDCType: {}", k)),
        },
        _ => Err(format!("missing LDCType: {}", k)),
    }
}



#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LDCType {
    Poly,
    Claret,
}

#[derive(Debug, Clone, Copy)]
pub struct LDC {
    ldc1: f64,
    ldc2: f64,
    ldc3: f64,
    ldc4: f64,
    mucrit: f64,
    ltype: LDCType,
}

impl LDC {
    // Default. Sets all to zero and type to POLY.
    pub fn new() -> Self {
        Self {
            ldc1: 0.0,
            ldc2: 0.0,
            ldc3: 0.0,
            ldc4: 0.0,
            mucrit: 0.0,
            ltype: LDCType::Poly,
        }
    }

    // Standard constructor
    pub fn with_params(
        ldc1: f64,
        ldc2: f64,
        ldc3: f64,
        ldc4: f64,
        mucrit: f64,
        ltype: LDCType,
    ) -> Self {
        Self {
            ldc1,
            ldc2,
            ldc3,
            ldc4,
            mucrit,
            ltype,
        }
    }

    /// Computes I(mu)
    pub fn imu(&self, mu: f64) -> f64 {
        if mu <= 0.0 {
            0.0
        } else {
            let mu = mu.min(1.0);
            let ommu = 1.0 - mu;
            let mut im = 1.0;

            match self.ltype {
                LDCType::Poly => {
                    im -= ommu
                        * (self.ldc1
                            + ommu * (self.ldc2 + ommu * (self.ldc3 + ommu * self.ldc4)));
                }
                LDCType::Claret => {
                    im -= self.ldc1 + self.ldc2 + self.ldc3 + self.ldc4;
                    let msq = mu.sqrt();
                    im += msq
                        * (self.ldc1
                            + msq * (self.ldc2 + msq * (self.ldc3 + msq * self.ldc4)));
                }
            }

            im
        }
    }

    /// To help applying mucrit
    pub fn see(&self, mu: f64) -> bool {
        mu > self.mucrit
    }
}


impl Default for LDC {
    fn default() -> Self {
        Self::new()
    }
}


pub type Etype = Vec<(f64, f64)>;


#[derive(Clone, Debug)]
pub struct Point {
    pub position: Vec3,
    pub direction: Vec3,
    pub area: f32,
    pub gravity: f32,
    pub eclipse: Etype,
    pub flux: f32,
}

impl Point {
    pub fn new(position: Vec3, direction: Vec3, area: f64, gravity: f64, eclipse: Etype) -> Self {
        Self {
            position,
            direction,
            area: area as f32,
            gravity: gravity as f32,
            eclipse,
            flux: 0.0,
        }
    }

    pub fn set_flux(&mut self, flux: f32) -> () {
        self.flux = flux;
    }

    pub fn is_visible(&self, phase: f64) -> bool {
        let phi: f64 = phase - phase.floor();
        for &(p1, p2) in &self.eclipse {
            if (phi >= p1 && phi <= p2) || phi <= p2 - 1.0 {
                return false
            }
        }
        true
    }
}

impl Default for Point {
    fn default() -> Self {
        Self::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            0.0,
            0.0,
            vec![(0.0, 0.0)],
            )
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Pparam {
    pub value: f64,
    pub range: f64,
    pub dstep: f64,
    pub vary: bool,
    pub defined: bool
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
        // let vary = fields.next().ok_or("missing vary")?.parse().map_err(|_| "bad vary")?;

        Ok(Pparam { value, range, dstep, vary: true, defined: true})
    }
}
#[derive(Debug, Clone, Copy)]
pub struct Model {

    pub q: Pparam,
    pub iangle: Pparam,
    pub r1: Pparam,
    pub r2: Pparam,
    pub cphi3: Pparam,
    pub cphi4: Pparam,
    pub spin1: Pparam,
    pub spin2: Pparam,
    pub t1: Pparam,
    pub t2: Pparam,
    pub ldc1_1: Pparam,
    pub ldc1_2: Pparam,
    pub ldc1_3: Pparam,
    pub ldc1_4: Pparam,
    pub ldc2_1: Pparam,
    pub ldc2_2: Pparam,
    pub ldc2_3: Pparam,
    pub ldc2_4: Pparam,
    pub velocity_scale: Pparam,
    pub beam_factor1: Pparam,
    pub beam_factor2: Pparam,
    pub t0: Pparam,
    pub period: Pparam,
    pub pdot: Pparam,
    pub deltat: Pparam,
    pub gravity_dark1: Pparam,
    pub gravity_dark2: Pparam,
    pub absorb: Pparam,
    pub slope: Pparam,
    pub quad: Pparam,
    pub cube: Pparam,
    pub third: Pparam,
    pub rdisc1: Pparam,
    pub rdisc2: Pparam,
    pub height_disc: Pparam,
    pub beta_disc: Pparam,
    pub temp_disc: Pparam,
    pub texp_disc: Pparam,
    pub lin_limb_disc: Pparam,
    pub quad_limb_disc: Pparam,
    pub temp_edge: Pparam,
    pub absorb_edge: Pparam,
    pub radius_spot: Pparam,
    pub length_spot: Pparam,
    pub height_spot: Pparam,
    pub expon_spot: Pparam,
    pub epow_spot: Pparam,
    pub angle_spot: Pparam,
    pub yaw_spot: Pparam,
    pub temp_spot: Pparam,
    pub tilt_spot: Pparam,
    pub cfrac_spot: Pparam,
    pub stsp11_long: Pparam,
    pub stsp11_lat: Pparam,
    pub stsp11_fwhm: Pparam,
    pub stsp11_tcen: Pparam,
    pub stsp12_long: Pparam,
    pub stsp12_lat: Pparam,
    pub stsp12_fwhm: Pparam,
    pub stsp12_tcen: Pparam,
    pub stsp13_long: Pparam,
    pub stsp13_lat: Pparam,
    pub stsp13_fwhm: Pparam,
    pub stsp13_tcen: Pparam,
    pub stsp21_long: Pparam,
    pub stsp21_lat: Pparam,
    pub stsp21_fwhm: Pparam,
    pub stsp21_tcen: Pparam,
    pub stsp22_long: Pparam,
    pub stsp22_lat: Pparam,
    pub stsp22_fwhm: Pparam,
    pub stsp22_tcen: Pparam,
    pub uesp_long1: Pparam,
    pub uesp_long2: Pparam,
    pub uesp_lathw: Pparam,
    pub uesp_taper: Pparam,
    pub uesp_temp: Pparam,
    pub delta_phase: f64,
    pub nlat1f: u32,
    pub nlat2f: u32,
    pub nlat1c: u32,
    pub nlat2c: u32,
    pub npole: bool,
    pub nlatfill: u32,
    pub nlngfill: u32,
    pub lfudge: f64,
    pub llo: f64,
    pub lhi: f64,
    pub phase1: f64,
    pub phase2: f64,
    pub wavelength: f64,
    pub roche1: bool,
    pub roche2: bool,
    pub eclipse1: bool,
    pub eclipse2: bool,
    pub glens1: bool,
    pub use_radii: bool,
    pub tperiod: f64,
    pub gdark_bolom1: bool,
    pub gdark_bolom2: bool,
    pub mucrit1: f64,
    pub mucrit2: f64,
    pub limb1: LDCType,
    pub limb2: LDCType,
    pub mirror: bool,
    pub add_disc: bool,
    pub nrad: u32,
    pub opaque: bool,
    pub add_spot: bool,
    pub nspot: u32,
    pub iscale: bool
}

impl Model {

    pub fn from_map(map: HashMap<String, Entry>) -> Result<Self, String> {
        Ok(Self {
            // Pparams
            q: get_p(&map, "q")?,
            iangle: get_p(&map, "iangle")?,
            r1: get_p(&map, "r1")?,
            r2: get_p(&map, "r2")?,
            cphi3: get_p(&map, "cphi3")?,
            cphi4: get_p(&map, "cphi4")?,
            spin1: get_p(&map, "spin1")?,
            spin2: get_p(&map, "spin2")?,
            t1: get_p(&map, "t1")?,
            t2: get_p(&map, "t2")?,
            ldc1_1: get_p(&map, "ldc1_1")?,
            ldc1_2: get_p(&map, "ldc1_2")?,
            ldc1_3: get_p(&map, "ldc1_3")?,
            ldc1_4: get_p(&map, "ldc1_4")?,
            ldc2_1: get_p(&map, "ldc2_1")?,
            ldc2_2: get_p(&map, "ldc2_2")?,
            ldc2_3: get_p(&map, "ldc2_3")?,
            ldc2_4: get_p(&map, "ldc2_4")?,
            velocity_scale: get_p(&map, "velocity_scale")?,
            beam_factor1: get_p(&map, "beam_factor1")?,
            beam_factor2: get_p(&map, "beam_factor2")?,
            t0: get_p(&map, "t0")?,
            period: get_p(&map, "period")?,
            pdot: get_p(&map, "pdot")?,
            deltat: get_p(&map, "deltat")?,
            gravity_dark1: get_p(&map, "gravity_dark1")?,
            gravity_dark2: get_p(&map, "gravity_dark2")?,
            absorb: get_p(&map, "absorb")?,
            slope: get_p(&map, "slope")?,
            quad: get_p(&map, "quad")?,
            cube: get_p(&map, "cube")?,
            third: get_p(&map, "third")?,
            rdisc1: get_p(&map, "rdisc1")?,
            rdisc2: get_p(&map, "rdisc2")?,
            height_disc: get_p(&map, "height_disc")?,
            beta_disc: get_p(&map, "beta_disc")?,
            temp_disc: get_p(&map, "temp_disc")?,
            texp_disc: get_p(&map, "texp_disc")?,
            lin_limb_disc: get_p(&map, "lin_limb_disc")?,
            quad_limb_disc: get_p(&map, "quad_limb_disc")?,
            temp_edge: get_p(&map, "temp_edge")?,
            absorb_edge: get_p(&map, "absorb_edge")?,
            radius_spot: get_p(&map, "radius_spot")?,
            length_spot: get_p(&map, "length_spot")?,
            height_spot: get_p(&map, "height_spot")?,
            expon_spot: get_p(&map, "expon_spot")?,
            epow_spot: get_p(&map, "epow_spot")?,
            angle_spot: get_p(&map, "angle_spot")?,
            yaw_spot: get_p(&map, "yaw_spot")?,
            temp_spot: get_p(&map, "temp_spot")?,
            tilt_spot: get_p(&map, "tilt_spot")?,
            cfrac_spot: get_p(&map, "cfrac_spot")?,

            stsp11_long: get_p(&map, "stsp11_long")?,
            stsp11_lat: get_p(&map, "stsp11_lat")?,
            stsp11_fwhm: get_p(&map, "stsp11_fwhm")?,
            stsp11_tcen: get_p(&map, "stsp11_tcen")?,

            stsp12_long: get_p(&map, "stsp12_long")?,
            stsp12_lat: get_p(&map, "stsp12_lat")?,
            stsp12_fwhm: get_p(&map, "stsp12_fwhm")?,
            stsp12_tcen: get_p(&map, "stsp12_tcen")?,

            stsp13_long: get_p(&map, "stsp13_long")?,
            stsp13_lat: get_p(&map, "stsp13_lat")?,
            stsp13_fwhm: get_p(&map, "stsp13_fwhm")?,
            stsp13_tcen: get_p(&map, "stsp13_tcen")?,

            stsp21_long: get_p(&map, "stsp21_long")?,
            stsp21_lat: get_p(&map, "stsp21_lat")?,
            stsp21_fwhm: get_p(&map, "stsp21_fwhm")?,
            stsp21_tcen: get_p(&map, "stsp21_tcen")?,

            stsp22_long: get_p(&map, "stsp22_long")?,
            stsp22_lat: get_p(&map, "stsp22_lat")?,
            stsp22_fwhm: get_p(&map, "stsp22_fwhm")?,
            stsp22_tcen: get_p(&map, "stsp22_tcen")?,

            uesp_long1: get_p(&map, "uesp_long1")?,
            uesp_long2: get_p(&map, "uesp_long2")?,
            uesp_lathw: get_p(&map, "uesp_lathw")?,
            uesp_taper: get_p(&map, "uesp_taper")?,
            uesp_temp: get_p(&map, "uesp_temp")?,
            
            // Scalars
            delta_phase: get_f64(&map, "delta_phase")?,
            nlat1f: get_u32(&map, "nlat1f")?,
            nlat2f: get_u32(&map, "nlat2f")?,
            nlat1c: get_u32(&map, "nlat1c")?,
            nlat2c: get_u32(&map, "nlat2c")?,
            npole: get_bool(&map, "npole")?,
            nlatfill: get_u32(&map, "nlatfill")?,
            nlngfill: get_u32(&map, "nlngfill")?,
            lfudge: get_f64(&map, "lfudge")?,
            llo: get_f64(&map, "llo")?,
            lhi: get_f64(&map, "lhi")?,
            phase1: get_f64(&map, "phase1")?,
            phase2: get_f64(&map, "phase2")?,
            wavelength: get_f64(&map, "wavelength")?,
            roche1: get_bool(&map, "roche1")?,
            roche2: get_bool(&map, "roche2")?,
            eclipse1: get_bool(&map, "eclipse1")?,
            eclipse2: get_bool(&map, "eclipse2")?,
            glens1: get_bool(&map, "glens1")?,
            use_radii: get_bool(&map, "use_radii")?,
            tperiod: get_f64(&map, "tperiod")?,
            gdark_bolom1: get_bool(&map, "gdark_bolom1")?,
            gdark_bolom2: get_bool(&map, "gdark_bolom2")?,
            mucrit1: get_f64(&map, "mucrit1")?,
            mucrit2: get_f64(&map, "mucrit2")?,
            limb1: get_ldc(&map, "limb1")?,
            limb2: get_ldc(&map, "limb2")?,
            mirror: get_bool(&map, "mirror")?,
            add_disc: get_bool(&map, "add_disc")?,
            nrad: get_u32(&map, "nrad")?,
            opaque: get_bool(&map, "opaque")?,
            add_spot: get_bool(&map, "add_spot")?,
            nspot: get_u32(&map, "nspot")?,
            iscale: get_bool(&map, "iscale")?,
        })
    }

    pub fn from_file(path: &str) -> Result<Self, String> {
        let map = load_entries(path)?;
        Self::from_map(map)
    }


    pub fn get_r1r2(&self) -> (f64, f64) {
        if self.use_radii {
            (self.r1.value, self.r2.value)
        }
        else {
            let sini = self.iangle.value.to_radians().sin();
            let r2pr1 = (1. - (sini * (2.*PI*self.cphi4.value).cos()).powi(2)).sqrt();
            let r2mr1 = (1. - (sini * (2.*PI*self.cphi3.value).cos()).powi(2)).sqrt();
            let rr1 = (r2pr1 - r2mr1)/2.;
            let rr2 = (r2pr1 + r2mr1)/2.;
            (rr1, rr2)
        }
    }

    pub fn get_ldc1(&self) -> LDC {
        LDC::with_params(self.ldc1_1.value, self.ldc1_2.value, self.ldc1_3.value, self.ldc1_4.value, self.mucrit1, self.limb1)
    }

    pub fn get_ldc2(&self) -> LDC {
        LDC::with_params(self.ldc2_1.value, self.ldc2_2.value, self.ldc2_3.value, self.ldc2_4.value, self.mucrit2, self.limb2)
    }
}