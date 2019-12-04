use super::fq2::Fq2;
use ff::{Field, PrimeField, PrimeFieldDecodingError, PrimeFieldRepr};

// B coefficient of BLS12-381 curve, 4.

pub const B_COEFF: Fq = Fq(FqRepr([
    0xcfdbf3f5f3fff9fb, 0xfe1fa0f000a06fff, 0xe8b7e7700e85fc3e, 0x44008ffadd3357e1, 0x282f3f0f67c4296c, 0x74c0425ecc910c4f, 0x502453a5c6d6e578, 0x0000000000001001,
]));

// The generators of G1/G2 are computed by finding the lexicographically smallest valid x coordinate,
// and its lexicographically smallest y coordinate and multiplying it by the cofactor such that the
// result is nonzero.

// Generator of G1
// x = 417884745340634798765233000978822987758734420362300207668066997064533934622737672342721452575663925077111609382014533036438905769338818242 
// y = 673637670223924330684198059886589534489270335351282802497505958433916523313083618690071972803140247298068633304810541127894772981811114826
pub const G1_GENERATOR_X: Fq = Fq(FqRepr([
    0x97a5db32331c2be4, 0x36c04af7dbbefe78, 0x1b8b0cca9ce8c0ca, 0x4d0393656da22a62, 0xc08fab50ecbbfd6d, 0x14d8ba84b3d406b9, 0xe9af7cbd7ba46b14, 0x0000000000000246
]));
pub const G1_GENERATOR_Y: Fq = Fq(FqRepr([
    0x1abfa1ee5e099134, 0x8a966c9845e9d8d3, 0x1ae92f5cacd181b5, 0x9b0fde0d3f6d22a9, 0xe642826e66b202b6, 0x3698de7949d17003, 0xf2343d75fec496f7, 0x0000000000000f94
]));

// Generator of G2
// x = u * 1846764664134969483257934161356810965175002412089114209667691682233401743148621727266658606012870874100477316617504704180537017742687710869 + 3707845102431748156503633645569191547125862856194755734098886781057541201234868851392057519860858682309588061502457035510910811039332431253
// y = u * 291882433578847694569900996612010244453092020179172559797602036034112117089451633355197146650399108661653511713316141585610054917105969330 + 2017449033840412018031501200847805342554627359533083916455379287540330689331615727730110078398984996391439558295252201121148377646992387980
pub const G2_GENERATOR_X_C0: Fq = Fq(FqRepr([
    0x0349c4f358b3b654, 0x7199173cae1ba450, 0x1251c85c811a8291, 0xa028d44c7ba987c9, 0x750d25b7434d63bd, 0xaa06e3e8054ef25d, 0x1cd83f6f0b5e2d07, 0x000000000000107a
]));
pub const G2_GENERATOR_X_C1: Fq = Fq(FqRepr([
    0x02117690d10230f6, 0xcbfbfc883c9f28de, 0x675332ecef7270f0, 0x995b99743435654b, 0x46391b2ab88c0e12, 0xbbf655794aca1fad, 0x328a296faa2861cb, 0x000000000000132b
]));
pub const G2_GENERATOR_Y_C0: Fq = Fq(FqRepr([
   0x3e65ea2fbec33f7e, 0xd269b22cea41b84f, 0x4411f6380e880fb6, 0x9dda5361e10039d1, 0x233df969e9bdf20b, 0x346df026c1170601, 0x4204a14dd9504209, 0x00000000000007cf
]));
pub const G2_GENERATOR_Y_C1: Fq = Fq(FqRepr([
    0x007ac198ce1bd854, 0x321a5d089b832efc, 0xccf88d4f2500bbf5, 0x23ef3b90f2828b1f, 0xb67263f562fc72c2, 0xa89643230e8bff96, 0xe2c768dacf428044, 0x0000000000000f4e
]));

// Coefficients for the Frobenius automorphism.
pub const FROBENIUS_COEFF: [[Fq2; 5];3] = [
        [Fq2 {c0: Fq(FqRepr([0x5ffb2514d3d275a9, 0x2baf303fc53f79ea, 0x8cc94daab46f576f, 0x437fc6ea62be8bea, 0xd9909df9285c1491, 0xb1a272261be8b53c, 0xff85964a34868db5, 0x0000000000000fee])),
        c1: Fq(FqRepr([0x4ab0daebd6d83502, 0xd470cfc59015d0c0, 0x825d5ce6a10fa8d0, 0x47c2e1f5979d90d5, 0x59481cd37c1efff3, 0x24f1a26d43d4b9e1, 0x55bfbf0325cecca0, 0x0000000000000566])),
        },
        Fq2 {c0: Fq(FqRepr([0x0, 0x0, 0x0, 0x0, 0x0, 0x0,0x0,0x0])),
        c1: Fq(FqRepr([0xed7c493b1f8e0750, 0x9511cd97fc844443, 0xe32591e007dbc947, 0x96ebac2c10cfb6e5, 0x510a96bc9f57be18, 0xa8c831720ceb8a5b, 0x3b14a327e8924c9c, 0x00000000000003e0])),
        },
        Fq2 {c0: Fq(FqRepr([0xc7a7e127e74e1fc5, 0xdb251c03bfa8f8d9, 0xabd17532a366f0c8, 0xd6becd6c10e79a86, 0xa8a5bb81331ba8bd, 0xf6855771a55a89d7, 0x83866a4a28890f53, 0x000000000000001e])),
        c1: Fq(FqRepr([0xc7a7e127e74e1fc5, 0xdb251c03bfa8f8d9, 0xabd17532a366f0c8, 0xd6becd6c10e79a86, 0xa8a5bb81331ba8bd, 0xf6855771a55a89d7, 0x83866a4a28890f53, 0x000000000000001e])),
        },
        Fq2 {c0: Fq(FqRepr([0x36c8463871e35b24, 0x5491b5d2a7570d99, 0x5989e117b61d8847, 0x851b25f2c98585ae, 0x8e6037cd502a0352, 0x50533ce4e82071a7, 0x39cc62be03b2af65, 0x000000000000028b])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        },
        Fq2 {c0: Fq(FqRepr([0x27a3063cbb20956e, 0x06d44c4384e872c4, 0x389ac2dd57d64838, 0x1a3e945673a62671, 0x8236597a5b77bd4f, 0xa827c997c1433f14, 0x830c00945d0f9d09, 0x000000000000100d])),
        c1: Fq(FqRepr([0x8308f9c3ef8a153d, 0xf94bb3c1d06cd7e6, 0xd68be7b3fda8b807, 0x7104148986b5f64e, 0xb0a2615249035735, 0x2e6c4afb9e7a3009, 0xd23954b8fd45bd4c, 0x0000000000000547])),
        }
        ],

        [Fq2 {c0: Fq(FqRepr([0xbd2fb6c58b1ca35b, 0x6b0e326d58d10666, 0x2c0118b14da336f8, 0xf456fcb3e98c65da, 0xe1ce24100523566b, 0x2dcbe32152d1e4c2, 0x1a30b22571c30db9, 0x0000000000001175])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        },
        Fq2 {c0: Fq(FqRepr([0x73e3b9c838c74f87, 0xab8e4a32adfe3d11, 0xb59cc9799f6177f8, 0x062782ed30d69711, 0xa47882ff54511132, 0x8640d7ae779cfd76, 0x1b78f28f56a2aaf0, 0x00000000000012ca])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        },
        Fq2 {c0: Fq(FqRepr([0xb6b40302adaaac2c, 0x408017c5552d36aa, 0x899bb0c851be4100, 0x11d08639474a3137,0xc2aa5eef4f2dbac6, 0x5874f48d24cb18b3, 0x01484069e4df9d37, 0x0000000000000155])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        },
        Fq2 {c0: Fq(FqRepr([0xed7c493b1f8e0750, 0x9511cd97fc844443, 0xe32591e007dbc947, 0x96ebac2c10cfb6e5, 0x510a96bc9f57be18, 0xa8c831720ceb8a5b, 0x3b14a327e8924c9c, 0x00000000000003e0])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        },
        Fq2 {c0: Fq(FqRepr([0x36c8463871e35b24, 0x5491b5d2a7570d99, 0x5989e117b61d8847, 0x851b25f2c98585ae, 0x8e6037cd502a0352, 0x50533ce4e82071a7, 0x39cc62be03b2af65, 0x000000000000028b])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        }
        ],

        [Fq2 {c0: Fq(FqRepr([0xe3041ed8c35c8ae6, 0x24fae40195ac51d0, 0x6355355eb2180f77, 0xb483db73e9748239, 0x8a32ff4b715f6bc6, 0xe00ebd21ba62e546, 0xd1beeb0331cc4b01, 0x0000000000001536])),
        c1: Fq(FqRepr([0xc7a7e127e74e1fc5, 0xdb251c03bfa8f8d9, 0xabd17532a366f0c8, 0xd6becd6c10e79a86, 0xa8a5bb81331ba8bd, 0xf6855771a55a89d7, 0x83866a4a28890f53, 0x000000000000001e])),
        },
        Fq2 {c0: Fq(FqRepr([0x0, 0x0, 0x0, 0x0, 0x0, 0x0,0x0, 0x0])),
        c1: Fq(FqRepr([0xf3f7fcfdfcfffe7f, 0xbf9fe840002813ff, 0x858af9c903c0bf3f, 0x797222a6b311eb88, 0x702e5bdd554d59be, 0x7e1f20063af2566a, 0x53fd14e37575bd1e, 0x0000000000001400])),
        },
        Fq2 {c0: Fq(FqRepr([0xe3041ed8c35c8ae6, 0x24fae40195ac51d0, 0x6355355eb2180f77, 0xb483db73e9748239, 0x8a32ff4b715f6bc6, 0xe00ebd21ba62e546, 0xd1beeb0331cc4b01, 0x0000000000001536])),
        c1: Fq(FqRepr([0xe3041ed8c35c8ae6, 0x24fae40195ac51d0, 0x6355355eb2180f77, 0xb483db73e9748239, 0x8a32ff4b715f6bc6, 0xe00ebd21ba62e546, 0xd1beeb0331cc4b01, 0x0000000000001536])),
        },
        Fq2 {c0: Fq(FqRepr([0xb6b40302adaaac2c, 0x408017c5552d36aa, 0x899bb0c851be4100, 0x11d08639474a3137, 0xc2aa5eef4f2dbac6, 0x5874f48d24cb18b3, 0x01484069e4df9d37, 0x0000000000000155])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0])),
        },
        Fq2 {c0: Fq(FqRepr([0xc7a7e127e74e1fc5, 0xdb251c03bfa8f8d9, 0xabd17532a366f0c8, 0xd6becd6c10e79a86, 0xa8a5bb81331ba8bd, 0xf6855771a55a89d7, 0x83866a4a28890f53, 0x000000000000001e])),
        c1: Fq(FqRepr([0xe3041ed8c35c8ae6, 0x24fae40195ac51d0, 0x6355355eb2180f77, 0xb483db73e9748239,0x8a32ff4b715f6bc6, 0xe00ebd21ba62e546, 0xd1beeb0331cc4b01,0x0000000000001536])),
        }
        ]
    ];

// -((2**384) mod q) mod q
pub const NEGATIVE_ONE: Fq = Fq(FqRepr([
    0xb6b40302adaaac2c, 0x408017c5552d36aa, 0x899bb0c851be4100, 0x11d08639474a3137, 0xc2aa5eef4f2dbac6, 0x5874f48d24cb18b3, 0x01484069e4df9d37, 0x0000000000000155
]));

// ((q - 3) / 4)
pub const P_3DIV4:[u64;8] = [
                0xaaab00002aaaaaaa, 0x00080001555552aa, 0x03c9aaa4555fc010, 0x22d0aa37fe970730, 0x8cb62eb3291ec521, 0x75a50524d7ef5bc7, 0x5551555356955695, 0x0000000000000555
            ];
// ((q - 1) / 2)
pub const P_1DIV2:[u64;8] = [
                0x5556000055555555, 0x00100002aaaaa555, 0x07935548aabf8020, 0x45a1546ffd2e0e60, 0x196c5d66523d8a42, 0xeb4a0a49afdeb78f, 0xaaa2aaa6ad2aad2a, 0x0000000000000aaa
            ];
                    
pub const COFACTORG1:[u64;3] = [0xaaa7fffeaaaaaaab, 0xffffd55aaab01556, 0x0000000001555554];

pub const COFACTORG2:[u64;10] = [0x8e371c70e38e38e5, 0x71e755538e31d553, 0x9d531bc31a9b7200, 0x130e799c48dc9183, 
                                 0xba89241e66e72c7a, 0x128ba285eba2329e, 0x574f227721f5f081, 0x41384ef449ef40f2, 0xf1d38e4555ca3436, 0x0000001c71c6ffff
                                ];

// The BLS parameter x for BLS12-461 is -0x1FFFFFFBFFFE00000000
pub const BLS_X_IS_NEGATIVE: bool = true;
// NAF Representation of the BLS x parameter end corresponding x div 2
pub const BLS_X_NAF: &'static [i8] = &[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

#[derive(PrimeField)]
#[PrimeFieldModulus = "3969508375500863470560772059146634051800057393085754326046523646985852496169198543994841284697713271737768244168253401239242781720740276907"]
#[PrimeFieldGenerator = "2"]
pub struct Fq(FqRepr);

 pub mod fr {
     use ff::{Field, PrimeField, PrimeFieldDecodingError, PrimeFieldRepr};
    #[derive(PrimeField)]
    #[PrimeFieldModulus = "521481194400158902870293791036394582812650143983424074083311820261824039635303638490268303361"]
    #[PrimeFieldGenerator = "2"]
    pub struct Fr(FrRepr);}
