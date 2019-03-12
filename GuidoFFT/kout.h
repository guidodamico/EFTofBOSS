#ifndef KOUT_H
#define KOUT_H

//const size_t Nout = 36 ; 
//const double kout[Nout] = { 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5 } ;

/*const size_t Nout = 80 ;
const double kout[Nout] = { 0.01, 0.01043818,0.01532369,0.02157276,0.02834865,0.03456784,0.04079704,
   0.0469901,0.05308914,0.05954472,0.06582334,0.07193966,0.07824595,
   0.08475764,0.09121436,0.0973611,0.1035747,0.1098226,0.116044,0.1223865,
   0.1286908,0.1349423,0.1411876,0.1475008,0.1537725,0.1600564,0.1663951,
   0.1727325,0.1790439,0.1853209,0.1915728,0.1978322,0.2041849,0.2104981,
   0.2167137,0.2229646,0.2292547,0.2355059,0.2417865,0.2480881,0.2543503,
   0.2606232,0.2669563,0.2732767,0.2796161,0.2858967,0.2921163,0.2983731,
   0.3046781,0.3109607,0.3172621,0.3235839,0.3298394,0.33614,0.3423984,
   0.3486832,0.3549478,0.3612031,0.3674905,0.3737731,0.3800737,0.3863721,
   0.3926236,0.3989204,0.4052014,0.4115043,0.4177986,0.4240402,0.4303546,
   0.4366253,0.4429172,0.4491835,0.4554587,0.4617794,0.4680951,0.4743754,
   0.4806578,0.4869284,0.4931946,0.4994603 } ;*/

const size_t Nout = 100 ;
const double kout[Nout] = {0.01, 0.01044427, 0.01090827, 0.01139289, 0.01189904,
       0.01242768, 0.0129798 , 0.01355645, 0.01415872, 0.01478775,
       0.01544472, 0.01613088, 0.01684752, 0.017596  , 0.01837773,
       0.0191942 , 0.02004693, 0.02093755, 0.02186774, 0.02283925,
       0.02385393, 0.02491368, 0.02602052, 0.02717652, 0.02838389,
       0.02964489, 0.03096192, 0.03233746, 0.03377411, 0.03527458,
       0.03684172, 0.03847848, 0.04018795, 0.04197337, 0.04383811,
       0.0457857 , 0.04781981, 0.04994429, 0.05216316, 0.0544806 ,
       0.05690099, 0.05942892, 0.06206916, 0.06482669, 0.06770673,
       0.07071472, 0.07385635, 0.07713755, 0.08056452, 0.08414374,
       0.08788197, 0.09178629, 0.09586406, 0.10012299, 0.10457113,
       0.10921689, 0.11406904, 0.11913676, 0.12442962, 0.12995763,
       0.13573122, 0.14176132, 0.14805932, 0.15, 0.15463712,
       0.15789474, 0.16150715, 0.16578947, 0.16868239, 0.17368421,
       0.17617641, 0.18157895, 0.18400335, 0.18947368, 0.19217803,
       0.19736842, 0.20071588, 0.20526316, 0.20963304, 0.21315789,
       0.21894636, 0.22105263, 0.22867344, 0.22894737, 0.23684211,
       0.23883266, 0.24473684, 0.24944323, 0.25263158, 0.26052518,
       0.26052632, 0.26842105, 0.27209948, 0.27631579, 0.28418798,
       0.28421053, 0.29210526, 0.29681353, 0.3, 0.31};


// const double kout[Nout] = { 0.01, 0.0104566, 0.010934, 0.0114332, 0.0119553, 0.0125011, 0.0130719, 0.0136688, 0.0142929, 0.0149454, 0.0156278, 0.0163414, 0.0170875, 0.0178677, 0.0186835, 0.0195366, 0.0204286, 0.0213613, 0.0223366, 0.0233565, 0.0244229, 0.025538, 0.0267041, 0.0279233, 0.0291983, 0.0305314, 0.0319254, 0.0333831, 0.0349073, 0.0365011, 0.0381677, 0.0399104, 0.0417327, 0.0436381, 0.0456306, 0.047714, 0.0498925, 0.0521706, 0.0545526, 0.0570434, 0.0596479, 0.0623713, 0.0652191, 0.0681969, 0.0713107, 0.0745666, 0.0779712, 0.0815313, 0.0852539, 0.0891464, 0.0932167, 0.0974729, 0.101923, 0.106577, 0.111443, 0.116531, 0.121852, 0.127416, 0.133233, 0.139317, 0.145678, 0.152329, 0.159284, 0.166557, 0.174162, 0.182113, 0.190428, 0.199123, 0.208215, 0.217722, 0.227662, 0.238057, 0.248927, 0.260292, 0.272177, 0.284604, 0.297598, 0.311186, 0.325395, 0.340252, 0.355787, 0.372032, 0.389018, 0.40678, 0.425353, 0.444774, 0.465082, 0.486317, 0.508521, 0.53174, 0.556018, 0.581405, 0.607951, 0.635709, 0.664735, 0.695086, 0.726822, 0.760008, 0.794709, 0.830994 } ;
/*const size_t Nout = 99 ;
const double kout[Nk] = { 5e-4, 1e-3, 4e-3, 7e-3 ,0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.67, 0.7, 0.73, 0.76, 0.79, 0.82, 0.85, 0.88, 0.91, 0.94, 0.97, 1., 1.04, 1.08, 1.12, 1.16, 1.2, 1.25, 1.3, 1.35, 1.4, 1.46, 1.52, 1.58, 1.65, 1.72, 1.79, 1.86, 1.94, 2.22, 2.3, 2.4, 2.5, 2.6, 2.8, 3., 3.2, 3.4, 3.7, 4., 4.4, 4.8, 5.3, 5.8, 6.4, 7., 7.5, 8., 8.6, 9.3, 10. } ;
*/

#endif
