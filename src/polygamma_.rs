const SQRHUGE: f64 = 1.3043817436596712000e+19;
const BRNTRM: usize = 9;
static BRN: [f64; 10] = [1.6666666666666666e-01, 3.3333333333333333e-02, 2.3809523809523809e-02, 3.3333333333333333e-02, 7.5757575757575757e-02, 2.5311355311355311e-01, 1.1666666666666667e+00, 7.0921568627450980e+00, 5.4971177944862155e+01, 5.2912424242424242e+02];

pub fn polygamma(k: usize, x: f64) -> Result<f64, &'static str>{
    if x < 0.0 {
        return Err("must be x > 0.0");
    }
    if x == 0.0 {
        return Ok(f64::INFINITY);
    }

    let mut s= 0.0;
    let mut slv = 13.06;

    if k > 3 {
        // calc slv
        let mut f = 1.0;
        for i in (21..(k+20)).rev() {
            f *= i as f64;
        }
        for i in (3..(k+2)).rev() {
            f /= i as f64;
        }
		f *= 174611.0 / 55.0;
		slv = 6.812921 * f.powf(1.0 / 18.0);
		if slv < 13.06 {
			slv = 13.06;
        }
    }

    let mut pk = 1.0;
    for i in 1..(k+1) {
        pk *= i as f64
    }

    if x >= slv {
        if x > SQRHUGE {
            return Ok(f64::INFINITY);
        }
        let x2 = x * x;
        let mut isgn = if k % 2 == 1 { -1 } else { 1 };
        if k == 0 {
            // digamma
            for i in (0..(BRNTRM+1)).rev() {
                let i2 = 2 * (i+1);
                s += BRN[i] / (i2 as f64) * (isgn as f64);
                s /= x2;
                isgn *= -1;
            }
            s += x.ln() - 0.5 / x;
        }else{
            // trigamma, tetragamma
            for i in (0..(BRNTRM+1)).rev() {
                let mut f = 1.0;
                let i2 = 2 * (i+1);
                let mut j = i2 + k - 1;
				while j > i2 {
                    f *= j as f64;
                    j -= 1;
                }
				s += BRN[i] * f * (isgn as f64);
				s /= x2;
                isgn *= -1; 
            }
            for _ in 0..k {
                s /= x; 
            }
            let mut pxk = 1.0;
            for _ in 0..k {
                pxk *= x;
            }

            s -= pk * 0.5 / pxk / x * (isgn as f64);
            let f = pk / (k as f64);
            s -= f / pxk * (isgn as f64);
        }
    }else{
		let n = (slv - x) as i64;
		let mut y = (n as f64) + x + 1.0;
		s = polygamma(k, y).unwrap();
		let isgn = if k % 2 == 1 { 1 } else { -1 };
        for _ in 0..(n+1) {
			y -= 1.0;
			if y.abs() < 1e-3 {
				if x > 0.0 {
					y = x - (((x + 0.5) as i64) as f64);
                }else {
					y = x - (((x - 0.5) as i64) as f64);
                }
			}
            let mut pxk = 1.0;
            for _ in 0..k {
                pxk *= y;
            }
			if pxk * y == 0.0 {
                return Ok(f64::INFINITY);
            }
			s += (isgn as f64) *  pk / pxk / y; 
		}
    }

    Ok(s)
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;
    const GAMMA: f64 = 0.577215664901532860606512090082402431042;
    const EPS: f64 = 1e-15;

    #[test]
    fn digamma_1_2() {
        assert!((polygamma(0, 0.5).unwrap() - (- GAMMA - 2.0 * (2.0f64).ln())).abs() < EPS);
    }

    #[test]
    fn digamma_1() {
        assert!((polygamma(0, 1.0).unwrap() - (- GAMMA )).abs() < EPS);
    }

    #[test]
    fn digamma_2() {
        assert!((polygamma(0, 2.0).unwrap() - (1.0 - GAMMA)).abs() < EPS);
    }

    #[test]
    fn digamma_3() {
        assert!((polygamma(0, 3.0).unwrap() - (- GAMMA + 1.5)).abs() < EPS);
    }

    #[test]
    fn digamma_4() {
        assert!((polygamma(0, 4.0).unwrap() - (- GAMMA + 11.0/6.0)).abs() < EPS);
    }

    #[test]
    fn trigamma_1_2() {
        assert!((polygamma(1, 0.5).unwrap() - (PI*PI/2.0)).abs() < EPS);
    }

    #[test]
    fn trigamma_1() {
        assert!((polygamma(1, 1.0).unwrap() - (PI*PI/6.0)).abs() < EPS);
    }

    #[test]
    fn trigamma_3_2() {
        assert!((polygamma(1, 1.5).unwrap() - (PI*PI/2.0 - 4.0)).abs() < EPS);
    }

    #[test]
    fn trigamma_2() {
        assert!((polygamma(1, 2.0).unwrap() - (PI*PI/6.0 - 1.0)).abs() < EPS);
    }

}