use ark_bn254::{fq::Fq, fq2::Fq2, G2Affine, G2Projective};
use ark_ff::{Field, BigInteger64};
use num_bigint::BigUint;
pub use sha2::{Sha256, digest::Digest};
use std::str::FromStr;
use crate::hash2g1;
use crate::hash2g1::Hash2FieldBN254;
use ark_ec::{AffineRepr, CurveGroup};

// MapToCurve2 implements the Shallue and van de Woestijne method, applicable to any elliptic curve in Weierstrass form
// No cofactor clearing or isogeny
// https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-16.html#straightline-svdw
#[cfg(all(feature = "constantine_compatible", not(feature = "gnark_crypto_compatible")))]
#[allow(non_snake_case)]
pub fn MapToCurve2(u: Fq2) -> G2Affine {

    // constants
    // https://github.com/mratsim/constantine/blob/master/constantine/named/constants/bn254_snarks_hash_to_curve_g2.nim

    let z = Fq2{
        c0: Fq::from_str("0").unwrap(), 
        c1: Fq::from_str("1").unwrap()
    };

    let c1 = Fq2{
        c0: Fq::from_str("19485874751759354771024239261021720505790618469301721065564631296452457478373").unwrap(), 
        c1: Fq::from_str("266929791119991161246907387137283842545076965332900288569378510910307636689").unwrap()
    };

    let c2 = Fq2{
        c0: Fq::from_str("0").unwrap() * Fq::from(Fq::R).inverse().unwrap(), 
        c1: Fq::from_str("10944121435919637611123202872628637544348155578648911831344518947322613104291").unwrap()
    };

    let c3 = Fq2{
        c0: Fq::from_str("8270257801618377462829664163334948115088143961679076698731296916415895764198").unwrap(), 
        c1: Fq::from_str("15403170217607925661891511707918230497750592932893890913125906786266381721360").unwrap()
    };

    let c4 = Fq2{
        c0: Fq::from_str("18685085378399381287283517099609868978155387573303020199856495763721534568303").unwrap(), 
        c1: Fq::from_str("355906388159988214995876516183045123393435953777200384759171347880410182252").unwrap()
    };

    let mut tv1 = u.square();       //    1.  tv1 = u²
    tv1 = tv1 * c1;                 //    2.  tv1 = tv1 * c1

    let tv2 = Fq2::ONE + tv1;       //    3.  tv2 = 1 + tv1

    tv1 = Fq2::ONE - tv1;           //    4.  tv1 = 1 - tv1
    let mut tv3 = tv1 * tv2;        //    5.  tv3 = tv1 * tv2

    tv3 = tv3.inverse().unwrap();   //    6.  tv3 = inv0(tv3)
    let mut tv4 = u * tv1;          //    7.  tv4 = u * tv1
    tv4 = tv4 * tv3;                //    8.  tv4 = tv4 * tv3
    tv4 = tv4 * c3;                 //    9.  tv4 = tv4 * c3
    let x1 = c2 - tv4;              //    10.  x1 = c2 - tv4

    let mut gx1 = x1.square();      //    11. gx1 = x1²
    //12. gx1 = gx1 + A     All curves in gnark-crypto have A=0 (j-invariant=0). It is crucial to include this step if the curve has nonzero A coefficient.
    gx1 = gx1 * x1;                 //    13. gx1 = gx1 * x1
    gx1 = gx1 + Fq2{c0: Fq::from_str("19485874751759354771024239261021720505790618469301721065564631296452457478373").unwrap(), c1: Fq::from_str("266929791119991161246907387137283842545076965332900288569378510910307636690").unwrap()};  //    14. gx1 = gx1 + B

    let x2 = c2 + tv4;              //    15.  x2 = c2 + tv4
    let mut gx2 = x2.square();      //    16. gx2 = x2²
    //    17. gx2 = gx2 + A (see 12.)
    gx2 = gx2 * x2;                 //    18. gx2 = gx2 * x2
    gx2 = gx2 + Fq2{c0: Fq::from_str("19485874751759354771024239261021720505790618469301721065564631296452457478373").unwrap(), c1: Fq::from_str("266929791119991161246907387137283842545076965332900288569378510910307636690").unwrap()};  //    19. gx2 = gx2 + B

    let mut x3 = tv2.square();      //    20.  x3 = tv2²
    x3 = x3 * tv3;                  //    21.  x3 = x3 * tv3
    x3 = x3.square();               //    22.  x3 = x3²
    x3 = x3 * c4;                   //    23.  x3 = x3 * c4

    x3 = x3 + z;                    //    24.  x3 = x3 + Z

    let mut x = if gx1.legendre().is_qr() {x1} else {x3};   //    25.   x = CMOV(x3, x1, e1)   # x = x1 if gx1 is square, else x = x3
    x = if gx2.legendre().is_qr() && !gx1.legendre().is_qr(){x2} else {x};      //    26.   x = CMOV(x, x2, e2)    # x = x2 if gx2 is square and gx1 is not

    let mut gx = x.square();        //    27.  gx = x²
    //    28.  gx = gx + A
    gx = gx * x;                    //    29.  gx = gx * x
    gx = gx + Fq2{c0: Fq::from_str("19485874751759354771024239261021720505790618469301721065564631296452457478373").unwrap(), c1: Fq::from_str("266929791119991161246907387137283842545076965332900288569378510910307636690").unwrap()};    //    30.  gx = gx + B

    let mut y = gx.sqrt().unwrap(); //    31.   y = sqrt(gx)

    #[allow(non_snake_case)]
    let signsNotEqual = g2Sgn0(u) ^ g2Sgn0(y);  //    32.  e3 = sgn0(u) == sgn0(y)
    tv1 = Fq2::ZERO - y;

    if signsNotEqual == 0 {y = y} else {y = tv1};   //    33.   y = CMOV(-y, y, e3)       # Select correct sign of y
    
    let res = G2Affine::new_unchecked(x, y);

    if !res.is_on_curve() {
        panic!("Point not on curve")
    }

    res

}

// MapToCurve2 implements the Shallue and van de Woestijne method, applicable to any elliptic curve in Weierstrass form
// No cofactor clearing or isogeny
// https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-16.html#straightline-svdw
#[cfg(all(feature = "gnark_crypto_compatible"))]
#[allow(non_snake_case)]
pub fn MapToCurve2(u: Fq2) -> G2Affine {

    //constants: https://github.com/Consensys/gnark-crypto/blob/master/ecc/bn254/hash_to_g2.go#L33
	//c1 = g(Z)
	//c2 = -Z / 2
	//c3 = sqrt(-g(Z) * (3 * Z² + 4 * A))     # sgn0(c3) MUST equal 0
	//c4 = -4 * g(Z) / (3 * Z² + 4 * A)

    let z = Fq2{
        c0: Fq::from_str("6350874878119819312338956282401532409788428879151445726012394534686998597021").unwrap() * Fq::from(Fq::R).inverse().unwrap(), 
        c1: Fq::from_str("0").unwrap()
    };

    let c1 = Fq2{
        c0: Fq::from_str("1234912246041461878588942434875861039904126177810565185887158306408069993214").unwrap() * Fq::from(Fq::R).inverse().unwrap(), 
        c1: Fq::from_str("568440292453150825972223760836185707764922522371208948902804025364325400423").unwrap() * Fq::from(Fq::R).inverse().unwrap()
    };

    let c2 = Fq2{
        c0: Fq::from_str("7768683996859727954953724731427871339453941139073188968338321679979113805781").unwrap() * Fq::from(Fq::R).inverse().unwrap(), 
        c1: Fq::from_str("0").unwrap()
    };

    let c3 = Fq2{
        c0: Fq::from_str("14301738157933195389348253840724448307870907218086206201704502222609770096511").unwrap() * Fq::from(Fq::R).inverse().unwrap(), 
        c1: Fq::from_str("18766576807588938823931941816656094168356257905513364070341807858523241306211").unwrap() * Fq::from(Fq::R).inverse().unwrap()
    };

    let c4 = Fq2{
        c0: Fq::from_str("12945612253170900976712347250337035339258705867784462193943147521219390814770").unwrap() * Fq::from(Fq::R).inverse().unwrap(), 
        c1: Fq::from_str("21130322481901740787616774064142360811676414460802878397485299194159459008019").unwrap() * Fq::from(Fq::R).inverse().unwrap()
    };

    let mut tv1 = u.square();       //    1.  tv1 = u²
    tv1 = tv1 * c1;                 //    2.  tv1 = tv1 * c1

    let tv2 = Fq2::ONE + tv1;       //    3.  tv2 = 1 + tv1

    tv1 = Fq2::ONE - tv1;           //    4.  tv1 = 1 - tv1
    let mut tv3 = tv1 * tv2;        //    5.  tv3 = tv1 * tv2

    tv3 = tv3.inverse().unwrap();   //    6.  tv3 = inv0(tv3)
    let mut tv4 = u * tv1;          //    7.  tv4 = u * tv1
    tv4 = tv4 * tv3;                //    8.  tv4 = tv4 * tv3
    tv4 = tv4 * c3;                 //    9.  tv4 = tv4 * c3
    let x1 = c2 - tv4;              //    10.  x1 = c2 - tv4

    let mut gx1 = x1.square();      //    11. gx1 = x1²
    //12. gx1 = gx1 + A     All curves in gnark-crypto have A=0 (j-invariant=0). It is crucial to include this step if the curve has nonzero A coefficient.
    gx1 = gx1 * x1;                 //    13. gx1 = gx1 * x1
    gx1 = gx1 + Fq2{c0: Fq::from_str("19485874751759354771024239261021720505790618469301721065564631296452457478373").unwrap(), c1: Fq::from_str("266929791119991161246907387137283842545076965332900288569378510910307636690").unwrap()};  //    14. gx1 = gx1 + B

    let x2 = c2 + tv4;              //    15.  x2 = c2 + tv4
    let mut gx2 = x2.square();      //    16. gx2 = x2²
    //    17. gx2 = gx2 + A (see 12.)
    gx2 = gx2 * x2;                 //    18. gx2 = gx2 * x2
    gx2 = gx2 + Fq2{c0: Fq::from_str("19485874751759354771024239261021720505790618469301721065564631296452457478373").unwrap(), c1: Fq::from_str("266929791119991161246907387137283842545076965332900288569378510910307636690").unwrap()};  //    19. gx2 = gx2 + B

    let mut x3 = tv2.square();      //    20.  x3 = tv2²
    x3 = x3 * tv3;                  //    21.  x3 = x3 * tv3
    x3 = x3.square();               //    22.  x3 = x3²
    x3 = x3 * c4;                   //    23.  x3 = x3 * c4

    x3 = x3 + z;                    //    24.  x3 = x3 + Z

    let mut x = if gx1.legendre().is_qr() {x1} else {x3};   //    25.   x = CMOV(x3, x1, e1)   # x = x1 if gx1 is square, else x = x3
    x = if gx2.legendre().is_qr() && !gx1.legendre().is_qr(){x2} else {x};      //    26.   x = CMOV(x, x2, e2)    # x = x2 if gx2 is square and gx1 is not

    let mut gx = x.square();        //    27.  gx = x²
    //    28.  gx = gx + A
    gx = gx * x;                    //    29.  gx = gx * x
    gx = gx + Fq2{c0: Fq::from_str("19485874751759354771024239261021720505790618469301721065564631296452457478373").unwrap(), c1: Fq::from_str("266929791119991161246907387137283842545076965332900288569378510910307636690").unwrap()};    //    30.  gx = gx + B

    let mut y = gx.sqrt().unwrap(); //    31.   y = sqrt(gx)

    #[allow(non_snake_case)]
    let signsNotEqual = g2Sgn0(u) ^ g2Sgn0(y);  //    32.  e3 = sgn0(u) == sgn0(y)
    tv1 = Fq2::ZERO - y;

    if signsNotEqual == 0 {y = y} else {y = tv1};   //    33.   y = CMOV(-y, y, e3)       # Select correct sign of y
    
    let res = G2Affine::new_unchecked(x, y);

    if !res.is_on_curve() {
        panic!("Point not on curve")
    }

    res

}

// g2Sgn0 is an algebraic substitute for the notion of sign in ordered fields
// Namely, every non-zero quadratic residue in a finite field of characteristic =/= 2 has exactly two square roots, one of each sign
// https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-16.html#name-the-sgn0-function
#[allow(non_snake_case)]
pub fn g2Sgn0(u: Fq2) -> u64 {
    let mut sign = 0u64;
    let mut zero = 1u64;

    let t: BigUint = u.c0.into();
    let mut sign_i = *BigUint::to_bytes_le(&t).get(0).unwrap() as u64 & 1;
    // let mut zero_i = hash2g1::g1NotZero(u.c0);
    let t: BigUint = u.c0.into();
    let zero_i = (*t.to_u64_digits().get(0).unwrap() == 0) as u64;
    sign = sign | (zero & sign_i);
    zero = zero & zero_i;

    let t: BigUint = u.c1.into();
    sign_i = *BigUint::to_bytes_le(&t).get(0).unwrap() as u64 & 1;

    sign = sign | (zero & sign_i);

    return sign;
}

#[allow(non_snake_case)]
pub fn g2NotZero(x: Fq2) -> u64 {
	//Assuming G1 is over Fp and that if hashing is available for G2, it also is for G1
	return hash2g1::g1NotZero(x.c0) | hash2g1::g1NotZero(x.c1);

}

// MapToG2 invokes the SVDW map, and guarantees that the result is in g2
#[allow(non_snake_case)]
pub fn MapToG2(u: Fq2) -> G2Affine {
	let res = MapToCurve2(u);
	ClearCofactor(res)
}

// HashToG2 hashes a message to a point on the G2 curve using the SVDW map.
// Slower than EncodeToG2, but usable as a random oracle.
// https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-16.html#roadmap
#[allow(non_snake_case)]
pub fn HashToG2(msg: &[u8], dst: &[u8]) -> G2Affine {
    let u = Fq::hash_to_field(msg, dst, 4);

    let q0 = MapToCurve2(
        Fq2{
            c0: u[0],
            c1: u[1],
        }
    );

    let q1 = MapToCurve2(
        Fq2{
            c0: u[2],
            c1: u[3],
        }
    );

    let q:G2Affine = (q0 + q1).into();

    ClearCofactor(q)
}

// https://github.com/Consensys/gnark-crypto/blob/master/ecc/bn254/g2.go#L635
#[allow(non_snake_case)]
pub fn ClearCofactor(q: G2Affine) -> G2Affine {
    
    const X_GEN: u64 = 4965661367192848881;

    let mut points = [G2Affine::identity();4];

    let x_gen_scalar = BigInteger64::from(X_GEN);

    points[0] = q.mul_bigint(x_gen_scalar).into();

    points[1] = (points[0] + points[0] + points[0]).into();

    points[1] = psi(&points[1]);

    points[2] = psi(&points[0]);
    points[2] = psi(&points[2]);

    points[3] = psi(&q);
    points[3] = psi(&points[3]);
    points[3] = psi(&points[3]);

    (points[0] + points[1] + points[2] + points[3]).into()

}


// ψ(p) = u o π o u⁻¹ where u:E'→E iso from the twist to E
pub fn psi(a: &G2Affine) -> G2Affine {
    
    let a: G2Projective = (*a).into();

    let mut p: G2Projective = G2Affine::identity().into();
    let endo_u = Fq2{
        c0: Fq::from_str("21575463638280843010398324269430826099269044274347216827212613867836435027261").unwrap(),
        c1: Fq::from_str("10307601595873709700152284273816112264069230130616436755625194854815875713954").unwrap()
    };

    let endo_v = Fq2{
        c0: Fq::from_str("2821565182194536844548159561693502659359617185244120367078079554186484126554").unwrap(),
        c1: Fq::from_str("3505843767911556378687030309984248845540243509899259641013678093033130930403").unwrap()
    };

    p.x = conjugate(&a.x);
    p.y = conjugate(&a.y);
    p.z = conjugate(&a.z);

    p.x = p.x * endo_u;
    p.y = p.y * endo_v;

    p.into_affine()
}


pub fn conjugate(a: &Fq2) -> Fq2 {

    let ax = a.c0;
    let ay = a.c1;

    Fq2::new(ax, -ay)
}

// EncodeToG2 hashes a message to a point on the G2 curve using the SVDW map.
// It is faster than HashToG2, but the result is not uniformly distributed. Unsuitable as a random oracle.
// https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-16.html#roadmap
#[allow(non_snake_case)]
pub fn EncodeToG2(msg: &[u8], dst: &[u8]) -> G2Affine {

    let u = Fq::hash_to_field(msg, dst, 2);
    let res = MapToCurve2(Fq2{
        c0: u[0],
        c1: u[1],
    });

    ClearCofactor(res)
}


// Test Vector: https://github.com/Consensys/gnark-crypto/blob/master/ecc/bn254/hash_vectors_test.go
#[cfg(all(feature = "gnark_crypto_compatible"))]
#[cfg(test)]
mod tests {

    use std::str::FromStr;
    use ark_bn254::Fq2;
    use ark_bn254::G2Affine;
    use crate::hash2g2::MapToCurve2;
    use crate::hash2g2::Fq;
    use crate::hash2g2::HashToG2;
    use crate::hash2g2::EncodeToG2;

    #[test]
    #[allow(non_snake_case)]
    fn MapToCurve2_test() {
        let u = Fq2{
            c0: Fq::from_str("15963713818282906360305918686195491545577210390832157279818305179904408824931").unwrap(),
            c1: Fq::from_str("2166278439352519416731010325104738631510195416620895094682522641528929475020").unwrap()
        };

        let q = MapToCurve2(u);
        let expected = G2Affine::new_unchecked(Fq2{
            c0: Fq::from_str("16872093352184426853297847012752141646605261411290781565485515569233955899058").unwrap(),
            c1: Fq::from_str("20482288690411193526247554560661659739533735966007371008469181348051437821826").unwrap()
        }, Fq2{
            c0: Fq::from_str("427035866446275812154335387235552457760650543923113579505536211797911740485").unwrap(),
            c1: Fq::from_str("14849552243024588631071292176876897701191437999604860450422231174965236442203").unwrap()
        });
        assert!(expected.is_on_curve());
        assert!(q == expected);

        let u = Fq2{
            c0: Fq::from_str("12752967732566665017975022503761080419696068755373050496264700974774108086129").unwrap(),
            c1: Fq::from_str("20655422394809824901799481664662586419100706577355794400212187554951433717414").unwrap()
        };

        let q = MapToCurve2(u);
        let expected = G2Affine::new_unchecked(Fq2{
            c0: Fq::from_str("12193882055337081757241417044229479753659926309860257758224177044622322698984").unwrap(),
            c1: Fq::from_str("10092155993942609715417531227866448864240630219985669320168414926220064901453").unwrap()
        }, Fq2{
            c0: Fq::from_str("21850450548984866542151665069165216760882062028063278212318726360439829725223").unwrap(),
            c1: Fq::from_str("10197523149668572844555341938160230574503097016636734560718180396672437043430").unwrap()
        });
        assert!(expected.is_on_curve());
        assert!(q == expected);

        let u = Fq2{
            c0: Fq::from_str("18898141882839095816276844526801422247849121311000147859768000750276893266433").unwrap(),
            c1: Fq::from_str("3788127287937052767604234353437582991385298973804519256517508390161626404924").unwrap()
        };

        let q = MapToCurve2(u);
        let expected = G2Affine::new_unchecked(Fq2{
            c0: Fq::from_str("452805888478466390914725495219599183584561454657558688011312346353060651482").unwrap(),
            c1: Fq::from_str("7959928416860499659800248632934402218020177178560427800377197797165640390130").unwrap()
        }, Fq2{
            c0: Fq::from_str("14268098188884406522254505541441598455366967966015814006726862271011081843493").unwrap(),
            c1: Fq::from_str("15148517265986515293057552799755027217326970615601185424102524485888012383276").unwrap()
        });
        assert!(expected.is_on_curve());
        assert!(q == expected);
    }


    #[test]
    #[allow(non_snake_case)]
    fn HashToCurve2_test() {
        let q = HashToG2(b"abc", b"QUUX-V01-CS02-with-BN254G2_XMD:SHA-256_SVDW_RO_");
        // println!("{:?}", q0);
        let expected = G2Affine::new_unchecked(Fq2{
            c0: Fq::from_str("10305213714312555419584685236164610766057227018997600762219755820581571775698").unwrap(),
            c1: Fq::from_str("5140998983273781645596043003996621170933075714207210952317183701750931672829").unwrap()
        }, Fq2{
            c0: Fq::from_str("12782657610222102886506935265351398708799194735435757564502179253917869011884").unwrap(),
            c1: Fq::from_str("15746452850775091549966312821847336261590899319279618339578671846526379873840").unwrap()
        });
        assert!(expected.is_on_curve());
        assert!(q.is_on_curve());
        assert!(q == expected);


        let q = HashToG2(b"", b"QUUX-V01-CS02-with-BN254G2_XMD:SHA-256_SVDW_RO_");
        let expected = G2Affine::new_unchecked(Fq2{
            c0: Fq::from_str("7947280525355502288245767042139433332619084425813891508679326584140902765312").unwrap(),
            c1: Fq::from_str("10530141512348869141982713319207053343182583313484148698392330696376288318261").unwrap()
        }, Fq2{
            c0: Fq::from_str("2079515028849057274649333561166551431956364880890028320215862191123161285080").unwrap(),
            c1: Fq::from_str("20169147323092870078028771345234445157617856249189458168875341276090072581620").unwrap()
        });
        assert!(expected.is_on_curve());
        assert!(q.is_on_curve());
        assert!(q == expected);

        let q = HashToG2(b"abcdef0123456789", b"QUUX-V01-CS02-with-BN254G2_XMD:SHA-256_SVDW_RO_");
        let expected = G2Affine::new_unchecked(Fq2{
            c0: Fq::from_str("9141649584568251133435811655082820452253999683001609355083509727807340928112").unwrap(),
            c1: Fq::from_str("19241337378620754008094815492162488101811979191715181531381201352430992486769").unwrap()
        }, Fq2{
            c0: Fq::from_str("18149222514336885092356998491550186845822771992585824025266466238465484336696").unwrap(),
            c1: Fq::from_str("9129360097802525322055823374454170177267012396640126715240529872313988489338").unwrap()
        });
        assert!(expected.is_on_curve());
        assert!(q.is_on_curve());
        assert!(q == expected);

        let q = HashToG2(b"q128_qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq", b"QUUX-V01-CS02-with-BN254G2_XMD:SHA-256_SVDW_RO_");
        let expected = G2Affine::new_unchecked(Fq2{
            c0: Fq::from_str("20353650816686918912609727598093385895712524005202794071238544969713808081729").unwrap(),
            c1: Fq::from_str("17684256473523682464984867199875609280081365245056171175421469718260504681254").unwrap()
        }, Fq2{
            c0: Fq::from_str("15896902550098660794387123920782326368527887924690142904247213645779094259076").unwrap(),
            c1: Fq::from_str("15390867031388969173331373188576779664345770454778413558467452103273727102977").unwrap()
        });
        assert!(expected.is_on_curve());
        assert!(q.is_on_curve());
        assert!(q == expected);

        let q = HashToG2(b"a512_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", b"QUUX-V01-CS02-with-BN254G2_XMD:SHA-256_SVDW_RO_");
        let expected = G2Affine::new_unchecked(Fq2{
            c0: Fq::from_str("16357539726107897952076989795377840344861047311782727672153303061989952217690").unwrap(),
            c1: Fq::from_str("10844839375884734385955874223756004111213539742547007380520745461640534925130").unwrap()
        }, Fq2{
            c0: Fq::from_str("20703414994053186684664027241143511234937261254193650036949701479117819278515").unwrap(),
            c1: Fq::from_str("11278285373922966720757356129051535273988981659843897570823718288010165493815").unwrap()
        });
        assert!(expected.is_on_curve());
        assert!(q.is_on_curve());
        assert!(q == expected);
    }

    #[test]
    fn encode_to_g2_test(){

        let q = EncodeToG2(b"abc", b"QUUX-V01-CS02-with-BN254G2_XMD:SHA-256_SVDW_NU_");
        let expected = G2Affine::new_unchecked(Fq2{
            c0: Fq::from_str("7290337032722028742894312496454770035215478865307401781131202361899492945880").unwrap(),
            c1: Fq::from_str("18605632812439984129247614998320701910992924251662446522071513278020164236983").unwrap()
        }, Fq2{
            c0: Fq::from_str("18565926515830203734257806009639634340842708214357641080049757818108383758101").unwrap(),
            c1: Fq::from_str("21026435153745179081072575128771379049563093023676092614267505429710510687357").unwrap()
        });
        assert!(expected.is_on_curve());
        assert!(q.is_on_curve());
        assert!(q == expected);

        let q = EncodeToG2(b"", b"QUUX-V01-CS02-with-BN254G2_XMD:SHA-256_SVDW_NU_");
        let expected = G2Affine::new_unchecked(Fq2{
            c0: Fq::from_str("2222545202255207121622252341720884612662004487208664408317925491033383016781").unwrap(),
            c1: Fq::from_str("3167015911722190124689644160541231412539898594125261078778351544051685395067").unwrap()
        }, Fq2{
            c0: Fq::from_str("20450065928984038040963910334909877834263207751235246619699259708122680403961").unwrap(),
            c1: Fq::from_str("4743914079645712786687283872900604142971897405422186449887746314931053675188").unwrap()
        });
        assert!(expected.is_on_curve());
        assert!(q.is_on_curve());
        assert!(q == expected);

        let q = EncodeToG2(b"abcdef0123456789", b"QUUX-V01-CS02-with-BN254G2_XMD:SHA-256_SVDW_NU_");
        let expected = G2Affine::new_unchecked(Fq2{
            c0: Fq::from_str("7148036967840401493869354348463445038937751410382870212181508408551260940454").unwrap(),
            c1: Fq::from_str("20374759774184409322905764368361574346849498692562411327726753719663647349306").unwrap()
        }, Fq2{
            c0: Fq::from_str("10483260641720359876745669935893009958901176103433938324656495809668720301952").unwrap(),
            c1: Fq::from_str("4967329811281913502786824686629199594924414673725274625361393684486574196665").unwrap()
        });
        assert!(expected.is_on_curve());
        assert!(q.is_on_curve());
        assert!(q == expected);

        let q = EncodeToG2(b"q128_qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq", b"QUUX-V01-CS02-with-BN254G2_XMD:SHA-256_SVDW_NU_");
        let expected = G2Affine::new_unchecked(Fq2{
            c0: Fq::from_str("13125957534682971537993382516266248139111688303006779518993454727694490448781").unwrap(),
            c1: Fq::from_str("17001778660286232066011321802406169530508086950847889879278742058506410997887").unwrap()
        }, Fq2{
            c0: Fq::from_str("21493575647416678969105094342755094397493417994062452257278523992404807814786").unwrap(),
            c1: Fq::from_str("17949331872970707337957877677444470806317301863003904868354967532347837920323").unwrap()
        });
        assert!(expected.is_on_curve());
        assert!(q.is_on_curve());
        assert!(q == expected);

        let q = EncodeToG2(b"a512_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", b"QUUX-V01-CS02-with-BN254G2_XMD:SHA-256_SVDW_NU_");
        let expected = G2Affine::new_unchecked(Fq2{
            c0: Fq::from_str("549777038834003445090108363667568008543458018879702424802185895608634070183").unwrap(),
            c1: Fq::from_str("17241878747514914600777537647804966909993790648100169250922892401106734291369").unwrap()
        }, Fq2{
            c0: Fq::from_str("8654939251469409238390702955149995024930650667376373847357704134187240104302").unwrap(),
            c1: Fq::from_str("3048816029845140188911702407210291715056845238089676496464503410870814591366").unwrap()
        });
        assert!(expected.is_on_curve());
        assert!(q.is_on_curve());
        assert!(q == expected);
    }
}

#[cfg(all(feature = "constantine_compatible", not(feature = "gnark_crypto_compatible")))]
#[cfg(test)]
mod tests {
    extern crate constantine_sys;
    
    use ark_bn254::G2Affine;
    use crate::hash2g2::HashToG2;
    use ark_bn254::G2Projective;
    use constantine_sys::*;
    use ark_ec::CurveGroup;
    use ::core::mem::MaybeUninit;
    use std::mem;
    // differential testing against constantine implementation: https://github.com/mratsim/constantine.git
    // https://github.com/mratsim/constantine/pull/437
    #[test]
    fn hash_to_curve_diff_test_g2(){

        // constantine output
        let mut result_constantine = MaybeUninit::<bn254_snarks_g2_jac>::uninit(); 
        let result_constantine_aff: G2Affine = unsafe {
            ctt_bn254_snarks_g2_jac_svdw_sha256(
                result_constantine.as_mut_ptr(),
                b"" as *const u8 ,
                0,
                b"abc" as *const u8,
                3,
                b"QUUX-V01-CS02-with-BN254G2_XMD:SHA-256_SVDW_RO_" as *const u8 ,
                47
            );
            
            let result_constantine_sys = mem::transmute::<MaybeUninit<bn254_snarks_g2_jac>, G2Projective>(result_constantine);
            result_constantine_sys.into_affine()
        };
        
        // native implementation output
        let result = HashToG2(b"abc", b"QUUX-V01-CS02-with-BN254G2_XMD:SHA-256_SVDW_RO_");

        //match the result
        assert_eq!(result, result_constantine_aff);
    }
}