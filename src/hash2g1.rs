
use ark_bn254::{fq::Fq, G1Affine as G1};
use ark_ff::Field;
use num_bigint::BigUint;
use digest::generic_array::GenericArray;
use num_integer::Integer;
use digest::generic_array::{typenum::U48, typenum::U32};
pub use sha2::{Sha256, digest::Digest};
use subtle::{Choice, ConditionallySelectable};
use std::str::FromStr;
pub trait From {
    fn from_bytes32(data: GenericArray::<u8, U32>) -> Self;
}

impl From for Fq {
    fn from_bytes32(bytes: GenericArray<u8, U32>) -> Self {
    Fq::from(BigUint::from_bytes_be(bytes.as_slice()))
    }
}

pub trait FromOkm<const L: usize>: Sized {
    /// Convert a byte sequence into a scalar
    fn from_okm(data: &[u8; L]) -> Self;
}

const L: usize = 48;
impl FromOkm<L> for Fq {
    fn from_okm(data: &[u8; L]) -> Self {
        let p = BigUint::from_bytes_be(
            &hex::decode("30644E72E131A029B85045B68181585D97816A916871CA8D3C208C16D87CFD47")
                .unwrap(),
        );

        let mut x = BigUint::from_bytes_be(&data[..]);
            x = x.mod_floor(&p);
            Fq::from(x)
    }
}


pub trait ExpandMsgSHA256 {
    /// Expands `msg` to the required number of bytes in `buf`
    #[allow(non_snake_case)]
    fn expand_message(msg: &[u8], dst: &[u8], LEN_IN_BYTES: usize) -> Vec<u8>;
}

#[allow(non_snake_case)]
impl ExpandMsgSHA256 for Fq {
    fn expand_message(msg: &[u8], dst: &[u8], LEN_IN_BYTES: usize) -> Vec<u8> {
        
        let b_in_bytes: usize = 32;
        let ell = (LEN_IN_BYTES + b_in_bytes - 1 )/ b_in_bytes;

        if ell > 255 {
            panic!("ell was too big in expand_message_xmd");
        }

        if dst.len() > 255 {
            panic!("dst size is invalid");
        }

        let b_0 = Sha256::new()
            .chain_update([0u8; 64])    // s_in_bytes for sha256 = 64
            .chain_update(msg)
            .chain_update([(LEN_IN_BYTES >> 8) as u8, LEN_IN_BYTES as u8, 0u8])
            .chain_update(dst)
            .chain_update([dst.len() as u8])
            .finalize();

        let mut b_vals = Sha256::new()
            .chain_update(&b_0[..])
            .chain_update([1u8])
            .chain_update(dst)
            .chain_update([dst.len() as u8])
            .finalize();

        let mut buf = [0u8; 4 * 48];
        let mut offset = 0;

        for i in 1..ell {
            // b_0 XOR b_(idx - 1)
            let mut tmp = GenericArray::<u8, U32>::default();
            b_0.iter()
                .zip(&b_vals[..])
                .enumerate()
                .for_each(|(j, (b0val, bi1val))| tmp[j] = b0val ^ bi1val);
            for b in b_vals {
                buf[offset % LEN_IN_BYTES].conditional_assign(
                    &b,
                    Choice::from(if offset < LEN_IN_BYTES { 1 } else { 0 }),
                );
                offset += 1;
            }
            b_vals = Sha256::new()
                .chain_update(tmp)
                .chain_update([(i + 1) as u8])
                .chain_update(dst)
                .chain_update([dst.len() as u8])
                .finalize();
        }
        for b in b_vals {
            buf[offset % LEN_IN_BYTES]
            .conditional_assign(&b, Choice::from(if offset < LEN_IN_BYTES { 1 } else { 0 }));
            offset += 1;
        }
        buf.into()
    }
}


pub trait Hash2FieldBN254 {
    fn hash_to_field (msg: &[u8], dst: &[u8], count: usize) -> Vec<Self> where Self: Sized;
}

impl Hash2FieldBN254 for Fq {
    
    fn hash_to_field(msg: &[u8], dst: &[u8], count: usize) -> Vec<Fq> {

        /*
        - p, the characteristic of F .
        - m, the extension degree of F, m >= 1 (see immediately above).
        - L = ceil((ceil(log2(p)) + k) / 8), where k is the security
            parameter of the suite (e.g., k = 128).
         */

        let len_per_elm = 48;
        let len_in_bytes = count * len_per_elm;
        // let len_in_bytes = count * len_per_elm;
        let pseudo_random_bytes = Fq::expand_message(msg, dst, len_in_bytes);
    
        let mut r = Vec::<Fq>::with_capacity(count);
        for i in 0..count {
            let bytes = GenericArray::<u8, U48>::from_slice(
                &pseudo_random_bytes[i * len_per_elm..(i + 1) * len_per_elm],
            );

            let x: [u8; 48] = bytes.as_slice().try_into().expect("Wrong length");
            r.push(Fq::from_okm(&x));
        }

        r
    }
}




// https://github.com/ConsenSys/gnark-crypto/blob/master/ecc/bn254/hash_to_g1.go
#[allow(non_snake_case)]
pub fn MapToCurve1(u: Fq) -> G1{

	//constants
	//c1 = g(Z)
	//c2 = -Z / 2
	//c3 = sqrt(-g(Z) * (3 * Z² + 4 * A))     # sgn0(c3) MUST equal 0
	//c4 = -4 * g(Z) / (3 * Z² + 4 * A)

	// Z := fp.Element{15230403791020821917, 754611498739239741, 7381016538464732716, 1011752739694698287}
	// c1 := fp.Element{1248766071674976557, 10548065924188627562, 16242874202584236114, 560012691975822483}
	// c2 := fp.Element{12997850613838968789, 14304628359724097447, 2950087706404981016, 1237622763554136189}
	// c3 := fp.Element{8972444824031832946, 5898165201680709844, 10690697896010808308, 824354360198587078}
	// c4 := fp.Element{12077013577332951089, 1872782865047492001, 13514471836495169457, 415649166299893576}

    let mut z: Fq = Fq::from_str("6350874878119819312338956282401532409788428879151445726012394534686998597021").unwrap();
    let mut c1: Fq = Fq::from_str("3515256640640002027109419384348854550457404359307959241360540244102768179501").unwrap();
    let mut c2: Fq = Fq::from_str("7768683996859727954953724731427871339453941139073188968338321679979113805781").unwrap();
    let mut c3: Fq = Fq::from_str("5174556184976127869173189452163337195348491024958816448391141365979064675186").unwrap();
    let mut c4: Fq = Fq::from_str("2609072103093089037936242735953952295622231240021995565748958972744717830193").unwrap();
    z = z * Fq::from(Fq::R).inverse().unwrap();
    c1 = c1 * Fq::from(Fq::R).inverse().unwrap();
    c2 = c2 * Fq::from(Fq::R).inverse().unwrap();
    c3 = c3 * Fq::from(Fq::R).inverse().unwrap();
    c4 = c4 * Fq::from(Fq::R).inverse().unwrap();

    let mut tv1: Fq = u.square();       //    1.  tv1 = u²
    tv1 = tv1 * c1;                     //    2.  tv1 = tv1 * c1
    let tv2: Fq = Fq::ONE + tv1;        //    3.  tv2 = 1 + tv1
    tv1 = Fq::ONE - tv1;                //    4.  tv1 = 1 - tv1
    let mut tv3: Fq = tv1 * tv2;        //    5.  tv3 = tv1 * tv2 
    
    tv3 = tv3.inverse().unwrap();       //    6.  tv3 = inv0(tv3)
    let mut tv4: Fq = u * tv1;          //    7.  tv4 = u * tv1  
    tv4 = tv4 * tv3;                    //    8.  tv4 = tv4 * tv3
    tv4 = tv4 * c3;                     //    9.  tv4 = tv4 * c3
    let x1: Fq = c2 - tv4;              //    10.  x1 = c2 - tv4
    
    let mut gx1: Fq = x1.square();      //    11. gx1 = x1²
    //12. gx1 = gx1 + A  It is crucial to include this step if the curve has nonzero A coefficient.
    gx1 = gx1 * x1;                     //    13. gx1 = gx1 * x1    
    gx1 = gx1 + Fq::from_str("3").unwrap();            //    14. gx1 = gx1 + B

    // let gx1NotSquare: i32 = if gx1.legendre().is_qr() {0} else {-1};    //    15.  e1 = is_square(gx1)
    // gx1NotSquare = 0 if gx1 is a square, -1 otherwise

    let x2: Fq = c2 + tv4;              //    16.  x2 = c2 + tv4
    let mut gx2: Fq = x2.square();      //    17. gx2 = x2²
    //    18. gx2 = gx2 + A     See line 12
    gx2 = gx2 * x2;                     //    19. gx2 = gx2 * x2
    gx2 = gx2 + Fq::from_str("3").unwrap();            //    20. gx2 = gx2 + B

    let mut x3: Fq = tv2.square();      //    22.  x3 = tv2²
    x3 = x3 * tv3;                      //    23.  x3 = x3 * tv3
    x3 = x3.square();                   //    24.  x3 = x3²
    x3 = x3 * c4;                       //    25.  x3 = x3 * c4

    x3 = x3 + z;                        //    26.  x3 = x3 + Z

    let mut x: Fq = if gx1.legendre().is_qr() { x1 } else { x3};  //    27.   x = CMOV(x3, x1, e1)   # x = x1 if gx1 is square, else x = x3

    if gx2.legendre().is_qr() && !gx1.legendre().is_qr() { x = x2 } //    28.   x = CMOV(x, x2, e2)    # x = x2 if gx2 is square and gx1 is not
	// Select x2 iff gx2 is square and gx1 is not, iff gx1SquareOrGx2Not = 0
    
    let mut gx = x.square();    //    29.  gx = x²
    //    30.  gx = gx + A
    gx = gx * x;                //    31.  gx = gx * x
    gx = gx + Fq::from_str("3").unwrap();      //    32.  gx = gx + B

    let mut y: Fq = gx.sqrt().unwrap();     //    33.   y = sqrt(gx)


    #[allow(non_snake_case)]
    let signsNotEqual = g1Sgn0(u) ^ g1Sgn0(y);

    tv1 = Fq::from(0) - y;
    //TODO: conditionallySelect
    if signsNotEqual == 0 {y = y} else {y = tv1}
    G1::new(x, y)


    
}

#[allow(non_snake_case)]
fn g1Sgn0(x: Fq) -> u64 {
    let t: BigUint = x.into();
    *BigUint::to_bytes_le(&t).get(0).unwrap() as u64 & 1
}

#[allow(non_snake_case)]
pub fn g1NotZero(x: Fq) -> u64 {

    let t: BigUint = x.into();
    let t_vec = t.to_u64_digits();
	return t_vec[0] | t_vec[1] | t_vec[2] | t_vec[3];

}

#[allow(non_snake_case)]
#[allow(dead_code)]
fn HashToG1(msg: &[u8], dst: &[u8]) -> G1 {
    let u = Hash2FieldBN254::hash_to_field(msg, dst, 2);
    let Q0 = MapToCurve1(u[0]);
    let Q1 = MapToCurve1(u[1]);
    let Q = Q0 + Q1;
    Q.into()
}

#[cfg(test)]
mod tests {

    use crate::hash2g1::Fq;
    use std::str::FromStr;
    use super::Hash2FieldBN254;
    use crate::hash2g1::MapToCurve1;
    use crate::hash2g1::G1;
    use crate::hash2g1::HashToG1;
    
    #[test]
    fn hash2field_test() {

        // Test Vector taken from https://github.com/Consensys/gnark-crypto/blob/master/ecc/bn254/hash_vectors_test.go

        // dst: []byte("QUUX-V01-CS02-with-BN254G1_XMD:SHA-256_SVDW_RO_"),
        // msg: "abc", P: point{"0x23f717bee89b1003957139f193e6be7da1df5f1374b26a4643b0378b5baf53d1", "0x4142f826b71ee574452dbc47e05bc3e1a647478403a7ba38b7b93948f4e151d"},
        // Q0: point{"0x1452c8cc24f8dedc25b24d89b87b64e25488191cecc78464fea84077dd156f8d", "0x209c3633505ba956f5ce4d974a868db972b8f1b69d63c218d360996bcec1ad41"},
        // Q1: point{"0x4e8357c98524e6208ae2b771e370f0c449e839003988c2e4ce1eaf8d632559f", "0x4396ec43dd8ec8f2b4a705090b5892219759da30154c39490fc4d59d51bb817"},
        // u0: "0x11945105b5e3d3b9392b5a2318409cbc28b7246aa47fa30da5739907737799a9", u1: "0x1255fc9ad5a6e0fb440916f091229bda611c41be2f2283c3d8f98c596be4c8c9",
        let u = Fq::hash_to_field(b"abc", b"QUUX-V01-CS02-with-BN254G1_XMD:SHA-256_SVDW_RO_", 2);
        assert!(u[0] == Fq::from_str("7951370986911800256774597109927097176311261202951929331835478768207980370345").unwrap());
        assert!(u[1] == Fq::from_str("8293556689416303717881563281438712057465092967957999993252567763605862533321").unwrap());

        // msg: "abcdef0123456789", P: point{"0x187dbf1c3c89aceceef254d6548d7163fdfa43084145f92c4c91c85c21442d4a", "0xabd99d5b0000910b56058f9cc3b0ab0a22d47cf27615f588924fac1e5c63b4d"},
        // Q0: point{"0x28d01790d2a1cc4832296774438acd46c2ce162d03099926478cf52319daba8d", "0x10227ab2707fd65fb45e87f0a48cfe3556f04113d27b1da9a7ae1709007355e1"},
        // Q1: point{"0x7dc256c7aadac1b4e1d23b3b2bbb5e2ffd9c753b9073d8d952ead8f812ce1b3", "0x2589008b2e15dcb3d16cdc1fed2634778001b1b28f0ab433f4f5ec6635c55e1e"},
        // u0: "0x2f7993a6b43a8dbb37060e790011a888157f456b895b925c3568690685f4983d", u1: "0x2677d0532b47a4cead2488845e7df7ebc16c0b8a2cd8a6b7f4ce99f51659794e",

        let u = Fq::hash_to_field(b"abcdef0123456789", b"QUUX-V01-CS02-with-BN254G1_XMD:SHA-256_SVDW_RO_", 2);
        assert!(u[0] == Fq::from_str("21473511429296129787161665655193361189518945362859158450118183976151186446397").unwrap());
        assert!(u[1] == Fq::from_str("17399580852346357386985693124899680967448413221719274165687915620563859110222").unwrap());

        // msg: "", P: point{"0xa976ab906170db1f9638d376514dbf8c42aef256a54bbd48521f20749e59e86", "0x2925ead66b9e68bfc309b014398640ab55f6619ab59bc1fab2210ad4c4d53d5"},
        // Q0: point{"0xe449b959abbd0e5ab4c873eaeb1ccd887f1d9ad6cd671fd72cb8d77fb651892", "0x29ff1e36867c60374695ee0c298fcbef2af16f8f97ed356fa75e61a797ebb265"},
        // Q1: point{"0x19388d9112a306fba595c3a8c63daa8f04205ad9581f7cf105c63c442d7c6511", "0x182da356478aa7776d1de8377a18b41e933036d0b71ab03f17114e4e673ad6e4"},
        // u0: "0x2f87b81d9d6ef05ad4d249737498cc27e1bd485dca804487844feb3c67c1a9b5", u1: "0x6de2d0d7c0d9c7a5a6c0b74675e7543f5b98186b5dbf831067449000b2b1f8e",

        let u = Fq::hash_to_field(b"", b"QUUX-V01-CS02-with-BN254G1_XMD:SHA-256_SVDW_RO_", 2);
        assert!(u[0] == Fq::from_str("21498498956904532351723378912032873852253513037650692457560050969314502748597").unwrap());
        assert!(u[1] == Fq::from_str("3106428082009635406807032300288584059640244342225966151234406580587112112014").unwrap());

        // msg: "a512_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", P: point{"0x1b05dc540bd79fd0fea4fbb07de08e94fc2e7bd171fe025c479dc212a2173ce", "0x1bf028afc00c0f843d113758968f580640541728cfc6d32ced9779aa613cd9b0"},
        // Q0: point{"0x2298ba379768da62495af6bb390ffca9156fde1dc167235b89c6dd008d2f2f3b", "0x660564cf6fce5cdea4780f5976dd0932559336fd072b4ddd83ec37f00fc7699"},
        // Q1: point{"0x2811dea430f7a1f6c8c941ecdf0e1e725b8ad1801ad15e832654bd8f10b62f16", "0x253390ed4fb39e58c30ca43892ab0428684cfb30b9df05fc239ab532eaa02444"},
        // u0: "0x48527470f534978bae262c0f3ba8380d7f560916af58af9ad7dcb6a4238e633", u1: "0x19a6d8be25702820b9b11eada2d42f425343889637a01ecd7672fbcf590d9ffe",

        let u = Fq::hash_to_field(b"a512_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", b"QUUX-V01-CS02-with-BN254G1_XMD:SHA-256_SVDW_RO_", 2);
        assert!(u[0] == Fq::from_str("2044513137826275527915612741016000753813717898656440700304636055936191489587").unwrap());
        assert!(u[1] == Fq::from_str("11602613730878338430727365363851039884306398846852682736694594518413917134846").unwrap());

    }

    #[test]
    #[allow(non_snake_case)]
    fn MapToCurve1_test() {
        let u = Fq::hash_to_field(b"abc", b"QUUX-V01-CS02-with-BN254G1_XMD:SHA-256_SVDW_RO_", 2);
        assert!(u[0] == Fq::from_str("7951370986911800256774597109927097176311261202951929331835478768207980370345").unwrap());
        assert!(u[1] == Fq::from_str("8293556689416303717881563281438712057465092967957999993252567763605862533321").unwrap());
        let q0 = MapToCurve1(u[0]);
        let q1 = MapToCurve1(u[1]);
        assert!(q0 == G1::new(Fq::from_str("9192524283969255398734814822241735402343760142215332184598869386265143635853").unwrap(), Fq::from_str("14750013374492649779039522357455217122947104756064249167130349093550158884161").unwrap()));
        assert!(q1 == G1::new(Fq::from_str("2219529064992744478098731193326567804904209297389738932911685687632211367327").unwrap(), Fq::from_str("1910726159786414357764375718946103460897900837832114831609513656424867805207").unwrap()));

        let u = Fq::hash_to_field(b"abcdef0123456789", b"QUUX-V01-CS02-with-BN254G1_XMD:SHA-256_SVDW_RO_", 2);
        assert!(u[0] == Fq::from_str("21473511429296129787161665655193361189518945362859158450118183976151186446397").unwrap());
        assert!(u[1] == Fq::from_str("17399580852346357386985693124899680967448413221719274165687915620563859110222").unwrap());
        let q0 = MapToCurve1(u[0]);
        let q1 = MapToCurve1(u[1]);
        assert!(q0 == G1::new(Fq::from_str("18460180777384996805517037410124907200489198402642233028065858702876325100173").unwrap(), Fq::from_str("7297925201307108404837100086863759533322513325723985709501528779399363778017").unwrap()));
        assert!(q1 == G1::new(Fq::from_str("3555154583542724794659651262588560064541528505277497563560719769602741821875").unwrap(), Fq::from_str("16977637197741440727690443467244845071598833410411827382713029829487302630942").unwrap()));

        let u = Fq::hash_to_field(b"", b"QUUX-V01-CS02-with-BN254G1_XMD:SHA-256_SVDW_RO_", 2);
        assert!(u[0] == Fq::from_str("21498498956904532351723378912032873852253513037650692457560050969314502748597").unwrap());
        assert!(u[1] == Fq::from_str("3106428082009635406807032300288584059640244342225966151234406580587112112014").unwrap());
        let q0 = MapToCurve1(u[0]);
        let q1 = MapToCurve1(u[1]);
        assert!(q0 == G1::new(Fq::from_str("6453599284581821454252898427469570073430843606970728650145294868078481709202").unwrap(), Fq::from_str("18995581315822946008285423533984677217009732542182181378734620089887646003813").unwrap()));
        assert!(q1 == G1::new(Fq::from_str("11407741707599100220112369632304941265828026024296299145123573579681208493329").unwrap(), Fq::from_str("10936143794657572576642578819087135925019845836839797797601194413922673415908").unwrap()));

        let u = Fq::hash_to_field(b"a512_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", b"QUUX-V01-CS02-with-BN254G1_XMD:SHA-256_SVDW_RO_", 2);
        assert!(u[0] == Fq::from_str("2044513137826275527915612741016000753813717898656440700304636055936191489587").unwrap());
        assert!(u[1] == Fq::from_str("11602613730878338430727365363851039884306398846852682736694594518413917134846").unwrap());
        let q0 = MapToCurve1(u[0]);
        let q1 = MapToCurve1(u[1]);
        assert!(q0 == G1::new(Fq::from_str("15648482829240231450830106706414350609765304380572340182587624553168656871227").unwrap(), Fq::from_str("2884090034870953736753092279678539049499634063513478774348194493913603274393").unwrap()));
        assert!(q1 == G1::new(Fq::from_str("18124086957708989159022363715051003299508967919164499453496316230748020813590").unwrap(), Fq::from_str("16826684847259413750643448296038340935081016395629901184485922600208379683908").unwrap()));

        // msg: "q128_qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq", P: point{"0xfe2b0743575324fc452d590d217390ad48e5a16cf051bee5c40a2eba233f5c", "0x794211e0cc72d3cbbdf8e4e5cd6e7d7e78d101ff94862caae8acbe63e9fdc78"},
        // Q0: point{"0x1c53b05f2fce15ba0b9100650c0fb46de1fb62f1d0968b69151151bd25dfefa4", "0x1fe783faf4bdbd79b717784dc59619106e4acccfe3b5d9750799729d855e7b81"},
        // Q1: point{"0x214a4e6e97adda47558f80088460eabd71ed35bc8ceafb99a493dd6f4e2b3f0a", "0xfaaeb29cc23f9d09b187a99741613aed84443e7c35736258f57982d336d13bd"},
        // u0: "0x2a50be15282ee276b76db1dab761f75401cdc8bd9fff81fcf4d428db16092a7b", u1: "0x23b41953676183c30aca54b5c8bd3ffe3535a6238c39f6b15487a5467d5d20eb",

        let u = Fq::hash_to_field(b"q128_qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq", b"QUUX-V01-CS02-with-BN254G1_XMD:SHA-256_SVDW_RO_", 2);
        assert!(u[0] == Fq::from_str("19139799307876008157674469077244497844490197231122854489816996874209678928507").unwrap());
        assert!(u[1] == Fq::from_str("16149156964295957170548772524136742336424608142546544142472739268994996707563").unwrap());
        let q0 = MapToCurve1(u[0]);
        let q1 = MapToCurve1(u[1]);
        assert!(q0 == G1::new(Fq::from_str("12812625340294489398993548898688649895606244530534785322892890577243153100708").unwrap(), Fq::from_str("14430750872577414903993343696062677117017041221676892932403418483169263778689").unwrap()));
        assert!(q1 == G1::new(Fq::from_str("15057612003824249181576746168110806738223995458659553230425471086211724164874").unwrap(), Fq::from_str("7086679767009137399570643369757025464320023320148085000688641996630730281917").unwrap()));

    }

    #[test]
    #[allow(non_snake_case)]
    fn HashToG1_test() {
        
        // Test Vector taken from https://github.com/Consensys/gnark-crypto/blob/master/ecc/bn254/hash_vectors_test.go

        let q = HashToG1(b"abc", b"QUUX-V01-CS02-with-BN254G1_XMD:SHA-256_SVDW_RO_");
        assert!(q == G1::new(Fq::from_str("16267524812466668166267883771992486438338357688076900798565538061554532963281").unwrap(), Fq::from_str("1844916233815282837483764409618609279507070495361570126601873459268232811805").unwrap()));

        let q = HashToG1(b"abcdef0123456789", b"QUUX-V01-CS02-with-BN254G1_XMD:SHA-256_SVDW_RO_");
        assert!(q == G1::new(Fq::from_str("11077683243901808951859264683654586764079462418577485658911541848692394044746").unwrap(), Fq::from_str("4858124309270455482359664916577923636817363175462672327824733704859450489677").unwrap()));

        let q = HashToG1(b"q128_qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq", b"QUUX-V01-CS02-with-BN254G1_XMD:SHA-256_SVDW_RO_");
        assert!(q == G1::new(Fq::from_str("449076125358095157945547407089359408531318284903480972761046551095956160348").unwrap(), Fq::from_str("3427911873443593747709927415036866402371639925174562008506349359915732032632").unwrap()));

        let q = HashToG1(b"", b"QUUX-V01-CS02-with-BN254G1_XMD:SHA-256_SVDW_RO_");
        assert!(q == G1::new(Fq::from_str("4790658965958450548702669593570794336562317867247372723806336874591549759110").unwrap(), Fq::from_str("1163238807669877429342450210709044731909255047583162173012265677391336920021").unwrap()));

        let q = HashToG1(b"a512_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", b"QUUX-V01-CS02-with-BN254G1_XMD:SHA-256_SVDW_RO_");
        assert!(q == G1::new(Fq::from_str("763925112321939766609678334678065587309331741428777416269918389033192485838").unwrap(), Fq::from_str("12636771015364464547273606234110225240317241569495907283228710706019336772016").unwrap()));

    }
}