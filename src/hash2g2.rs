use ark_bn254::{fq::Fq, fq2::Fq2, G2Affine, G1Projective, G2Projective};
// use ark_ec::short_weierstrass::Affine;
use ark_ff::Field;
use num_bigint::BigUint;
// use digest::generic_array::GenericArray;
// use num_integer::Integer;
// use digest::generic_array::{typenum::U48, typenum::U32};
pub use sha2::{Sha256, digest::Digest};
// use subtle::{Choice, ConditionallySelectable};
use std::str::FromStr;
use crate::hash2g1;
// use crate::hash2g1::ExpandMsgSHA256;
use crate::hash2g1::Hash2FieldBN254;
// use crate::hash2g1::FromOkm;
use ark_ec::AffineRepr;

#[allow(non_snake_case)]
pub fn MapToCurve2(u: Fq2) -> G2Affine {

    //constants
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

    let mut tv1 = u.square();
    tv1 = tv1 * c1;

    let tv2 = Fq2::ONE + tv1;

    tv1 = Fq2::ONE - tv1;
    let mut tv3 = tv1 * tv2;

    tv3 = tv3.inverse().unwrap();
    let mut tv4 = u * tv1;
    tv4 = tv4 * tv3;
    tv4 = tv4 * c3;
    let x1 = c2 - tv4;

    let mut gx1 = x1.square();
    gx1 = gx1 * x1;
    gx1 = gx1 + Fq2{c0: Fq::from_str("19485874751759354771024239261021720505790618469301721065564631296452457478373").unwrap(), c1: Fq::from_str("266929791119991161246907387137283842545076965332900288569378510910307636690").unwrap()};

    let x2 = c2 + tv4;
    let mut gx2 = x2.square();
    gx2 = gx2 * x2;
    gx2 = gx2 + Fq2{c0: Fq::from_str("19485874751759354771024239261021720505790618469301721065564631296452457478373").unwrap(), c1: Fq::from_str("266929791119991161246907387137283842545076965332900288569378510910307636690").unwrap()};

    let mut x3 = tv2.square();
    x3 = x3 * tv3;
    x3 = x3.square();
    x3 = x3 * c4;

    x3 = x3 + z;

    let mut x = if gx1.legendre().is_qr() {x1} else {x3};
    x = if gx2.legendre().is_qr() && !gx1.legendre().is_qr(){x2} else {x};

    let mut gx = x.square();
    gx = gx * x;
    gx = gx + Fq2{c0: Fq::from_str("19485874751759354771024239261021720505790618469301721065564631296452457478373").unwrap(), c1: Fq::from_str("266929791119991161246907387137283842545076965332900288569378510910307636690").unwrap()};

    let mut y = gx.sqrt().unwrap();

    #[allow(non_snake_case)]
    let signsNotEqual = g2Sgn0(u) ^ g2Sgn0(y);
    tv1 = Fq2::ZERO - y;

    if signsNotEqual == 0 {y = y} else {y = tv1};
    
    let res = G2Affine::new_unchecked(x, y);

    if !res.is_on_curve() {
        panic!("Point not on curve")
    }

    res

}

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


#[allow(non_snake_case)]
pub fn HashToG2(msg: &[u8], dst: &[u8]) -> G2Affine {
    let u = Fq::hash_to_field(msg, dst, 4);

    let q0 = MapToCurve2(
        Fq2{
            c0: u[0],
            c1: u[1],
        }
    );

    // println!("{:?}", q0.x.c0 * Fq::from(Fq::R));

    let q1 = MapToCurve2(
        Fq2{
            c0: u[2],
            c1: u[3],
        }
    );


    let q:G2Affine = (q0 + q1).into();
    // println!("{:?}", q.y.c0 * Fq::from(Fq::R));

    let _q = q.clear_cofactor();
    // let qr: G2Projective = _q.into();
    println!("{:?}", _q.x.c0 );
    _q.into()
}

#[cfg(test)]
mod tests {

    use std::str::FromStr;
    use ark_bn254::Fq2;
    use ark_bn254::G2Affine;
    use crate::hash2g2::MapToCurve2;
    use crate::hash2g2::Fq;
    use crate::hash2g2::HashToG2;
    // use ark_ff::Field;
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


    // #[test]
    // #[allow(non_snake_case)]
    // fn HashToCurve2_test() {

    //     let q0 = HashToG2(b"abc", b"QUUX-V01-CS02-with-BN254G2_XMD:SHA-256_SVDW_RO_");
    //     // println!("{:?}", q0);
    //     let expected = G2Affine::new_unchecked(Fq2{
    //         c0: Fq::from_str("10305213714312555419584685236164610766057227018997600762219755820581571775698").unwrap(),
    //         c1: Fq::from_str("5140998983273781645596043003996621170933075714207210952317183701750931672829").unwrap()
    //     }, Fq2{
    //         c0: Fq::from_str("12782657610222102886506935265351398708799194735435757564502179253917869011884").unwrap(),
    //         c1: Fq::from_str("15746452850775091549966312821847336261590899319279618339578671846526379873840").unwrap()
    //     });
    //     assert!(expected.is_on_curve());
    //     assert!(q0.x.c0 == expected.x.c0);

        // let u = Fq2{
        //     c0: Fq::from_str("12752967732566665017975022503761080419696068755373050496264700974774108086129").unwrap(),
        //     c1: Fq::from_str("20655422394809824901799481664662586419100706577355794400212187554951433717414").unwrap()
        // };

        // let q = MapToCurve2(u);
        // let expected = G2Affine::new_unchecked(Fq2{
        //     c0: Fq::from_str("12193882055337081757241417044229479753659926309860257758224177044622322698984").unwrap(),
        //     c1: Fq::from_str("10092155993942609715417531227866448864240630219985669320168414926220064901453").unwrap()
        // }, Fq2{
        //     c0: Fq::from_str("21850450548984866542151665069165216760882062028063278212318726360439829725223").unwrap(),
        //     c1: Fq::from_str("10197523149668572844555341938160230574503097016636734560718180396672437043430").unwrap()
        // });
        // assert!(expected.is_on_curve());
        // assert!(q.x == expected.x);

        // let u = Fq2{
        //     c0: Fq::from_str("18898141882839095816276844526801422247849121311000147859768000750276893266433").unwrap(),
        //     c1: Fq::from_str("3788127287937052767604234353437582991385298973804519256517508390161626404924").unwrap()
        // };

        // let q = MapToCurve2(u);
        // let expected = G2Affine::new_unchecked(Fq2{
        //     c0: Fq::from_str("452805888478466390914725495219599183584561454657558688011312346353060651482").unwrap(),
        //     c1: Fq::from_str("7959928416860499659800248632934402218020177178560427800377197797165640390130").unwrap()
        // }, Fq2{
        //     c0: Fq::from_str("14268098188884406522254505541441598455366967966015814006726862271011081843493").unwrap(),
        //     c1: Fq::from_str("15148517265986515293057552799755027217326970615601185424102524485888012383276").unwrap()
        // });
        // assert!(expected.is_on_curve());
        // assert!(q.x == expected.x);
    // }

}