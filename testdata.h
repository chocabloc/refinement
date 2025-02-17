#include <bits/stdc++.h>
#include "fgdata.h"
#include "Eigen/Core"

static Eigen::VectorXd test_rawX {{
    -0.396590, 3.112523, 0.042769,
    -0.412446, 3.799885, -0.162635,
    0.394289, 1.949852, 0.322617,
    0.790872, 2.012496, 0.252919,
    1.555705, 1.887771, 0.460957,
    2.207440, 2.689131, 0.177323,
    1.163043, 1.940795, 0.615613,
    1.685810, 2.315545, 0.396239,
    2.375365, 0.685177, 0.772358,
    2.427519, 0.019110, 1.957537,
    3.262224, 0.008293, -0.114206,
    3.638820, 0.198608, -1.442452,
    3.817412, -1.060165, 0.599924,
    3.293521, -1.031139, 1.861131,
    1.867268, 0.291299, 2.834295,
    3.236563, 1.005719, -2.024873,
    4.544349, -0.666524, -2.007220,
    4.729306, -1.931262, 0.029251,
    3.499935, -1.658956, 2.578785,
    4.847888, -0.534550, -3.032427,
    5.080429, -1.721333, -1.272739,
    5.150675, -2.749541, 0.581926,
    5.784621, -2.372144, -1.754142,
    -0.467663, 0.721431, 0.569207,
    -0.226113, -0.069734, 1.444112,
    -1.473785, 0.568046, -0.207598,
    -1.626186, 1.225359, -0.891086,
    -2.373377, -0.563858, -0.075087,
    -2.023278, -1.261404, -0.503222,
    -3.790972, -0.156450, -0.461124,
    -3.936923, 0.380162, -0.727657,
    -3.935045, 0.379162, -0.726935,
    -3.935352, 0.379870, -0.727303,
    -4.823357, -1.204617, -0.125684,
    -4.492305, -2.246835, 0.398670,
    -6.082520, -0.932044, -0.427635,
    -6.272495, -0.081935, -0.844084,
    -6.770367, -1.592341, -0.221593,
    -6.521425, -0.837166, -0.532825,
    -2.379758, -1.210041, 1.013368,
    -3.413096, -1.688194, 1.373527,
}};
static DistConstraints test_equality_cons = {
    DistConstraint(3, 1, 1, 1.458563),
    DistConstraint(24, 1, 1, 2.457634),
    DistConstraint(1, 2, 1, 0.980362),
    DistConstraint(3, 2, 1, 2.117337),
    DistConstraint(1, 3, 1, 1.458563),
    DistConstraint(2, 3, 1, 2.117337),
    DistConstraint(4, 3, 1, 1.080082),
    DistConstraint(5, 3, 1, 1.529548),
    DistConstraint(6, 3, 1, 2.141912),
    DistConstraint(24, 3, 1, 1.524672),
    DistConstraint(26, 3, 1, 2.424892),
    DistConstraint(1, 4, 1, 2.063538),
    DistConstraint(2, 4, 1, 2.275314),
    DistConstraint(3, 4, 1, 1.080082),
    DistConstraint(24, 4, 1, 2.133947),
    DistConstraint(1, 5, 1, 2.452413),
    DistConstraint(3, 5, 1, 1.529548),
    DistConstraint(4, 5, 1, 2.140024),
    DistConstraint(6, 5, 1, 1.079662),
    DistConstraint(14, 5, 1, 3.698906),
    DistConstraint(16, 5, 1, 3.138995),
    DistConstraint(20, 5, 1, 5.401299),
    DistConstraint(3, 6, 1, 2.141912),
    DistConstraint(5, 6, 1, 1.079662),
    DistConstraint(3, 7, 1, 2.141270),
    DistConstraint(5, 7, 1, 1.080490),
    DistConstraint(6, 7, 1, 1.761470),
    DistConstraint(3, 8, 1, 1.952422),
    DistConstraint(5, 8, 1, 0.625583),
    DistConstraint(6, 8, 1, 0.879951),
    DistConstraint(7, 8, 1, 0.881519),
    DistConstraint(3, 9, 1, 2.530408),
    DistConstraint(5, 9, 1, 1.498787),
    DistConstraint(6, 9, 1, 2.100228),
    DistConstraint(7, 9, 1, 2.101292),
    DistConstraint(8, 9, 1, 1.907492),
    DistConstraint(5, 10, 1, 2.562847),
    DistConstraint(9, 10, 1, 1.365692),
    DistConstraint(5, 11, 1, 2.620101),
    DistConstraint(9, 11, 1, 1.431728),
    DistConstraint(10, 11, 1, 2.239708),
    DistConstraint(5, 12, 1, 3.305122),
    DistConstraint(9, 12, 1, 2.604164),
    DistConstraint(10, 12, 1, 3.623475),
    DistConstraint(11, 12, 1, 1.397395),
    DistConstraint(5, 13, 1, 3.744450),
    DistConstraint(9, 13, 1, 2.286267),
    DistConstraint(10, 13, 1, 2.233189),
    DistConstraint(11, 13, 1, 1.408670),
    DistConstraint(12, 13, 1, 2.415417),
    DistConstraint(5, 14, 1, 3.698906),
    DistConstraint(9, 14, 1, 2.244375),
    DistConstraint(10, 14, 1, 1.374086),
    DistConstraint(11, 14, 1, 2.240519),
    DistConstraint(12, 14, 1, 3.552953),
    DistConstraint(13, 14, 1, 1.369773),
    DistConstraint(20, 14, 1, 5.172214),
    DistConstraint(5, 15, 1, 2.889417),
    DistConstraint(9, 15, 1, 2.165711),
    DistConstraint(10, 15, 1, 1.079232),
    DistConstraint(11, 15, 1, 3.283564),
    DistConstraint(12, 15, 1, 4.642863),
    DistConstraint(13, 15, 1, 3.273246),
    DistConstraint(14, 15, 1, 2.187248),
    DistConstraint(5, 16, 1, 3.138995),
    DistConstraint(9, 16, 1, 2.952158),
    DistConstraint(10, 16, 1, 4.193533),
    DistConstraint(11, 16, 1, 2.163353),
    DistConstraint(12, 16, 1, 1.080065),
    DistConstraint(13, 16, 1, 3.406174),
    DistConstraint(14, 16, 1, 4.403864),
    DistConstraint(15, 16, 1, 5.112386),
    DistConstraint(20, 16, 1, 2.460538),
    DistConstraint(5, 17, 1, 4.666398),
    DistConstraint(9, 17, 1, 3.791009),
    DistConstraint(10, 17, 1, 4.560445),
    DistConstraint(11, 17, 1, 2.392452),
    DistConstraint(12, 17, 1, 1.381863),
    DistConstraint(13, 17, 1, 2.742408),
    DistConstraint(14, 17, 1, 4.092765),
    DistConstraint(15, 17, 1, 5.632313),
    DistConstraint(16, 17, 1, 2.137909),
    DistConstraint(5, 18, 1, 5.018720),
    DistConstraint(9, 18, 1, 3.620890),
    DistConstraint(10, 18, 1, 3.599205),
    DistConstraint(11, 18, 1, 2.453378),
    DistConstraint(12, 18, 1, 2.826496),
    DistConstraint(13, 18, 1, 1.392306),
    DistConstraint(14, 18, 1, 2.505504),
    DistConstraint(15, 18, 1, 4.604443),
    DistConstraint(16, 18, 1, 3.906341),
    DistConstraint(17, 18, 1, 2.414061),
    DistConstraint(5, 19, 1, 4.594912),
    DistConstraint(9, 19, 1, 3.184722),
    DistConstraint(10, 19, 1, 2.100508),
    DistConstraint(11, 19, 1, 3.189016),
    DistConstraint(12, 19, 1, 4.446883),
    DistConstraint(13, 19, 1, 2.097817),
    DistConstraint(14, 19, 1, 0.980431),
    DistConstraint(15, 19, 1, 2.573866),
    DistConstraint(16, 19, 1, 5.346336),
    DistConstraint(17, 19, 1, 4.820218),
    DistConstraint(18, 19, 1, 2.851811),
    DistConstraint(5, 20, 1, 5.401299),
    DistConstraint(9, 20, 1, 4.715100),
    DistConstraint(10, 20, 1, 5.589857),
    DistConstraint(11, 20, 1, 3.375743),
    DistConstraint(12, 20, 1, 2.136107),
    DistConstraint(13, 20, 1, 3.822320),
    DistConstraint(14, 20, 1, 5.172214),
    DistConstraint(15, 20, 1, 6.651959),
    DistConstraint(16, 20, 1, 2.460538),
    DistConstraint(17, 20, 1, 1.080193),
    DistConstraint(18, 20, 1, 3.378752),
    DistConstraint(19, 20, 1, 5.895418),
    DistConstraint(5, 21, 1, 5.367550),
    DistConstraint(9, 21, 1, 4.180972),
    DistConstraint(10, 21, 1, 4.546663),
    DistConstraint(11, 21, 1, 2.780117),
    DistConstraint(12, 21, 1, 2.423946),
    DistConstraint(13, 21, 1, 2.362030),
    DistConstraint(14, 21, 1, 3.684728),
    DistConstraint(15, 21, 1, 5.611944),
    DistConstraint(16, 21, 1, 3.400430),
    DistConstraint(17, 21, 1, 1.401235),
    DistConstraint(18, 21, 1, 1.368423),
    DistConstraint(19, 21, 1, 4.175128),
    DistConstraint(20, 21, 1, 2.144241),
    DistConstraint(5, 22, 1, 5.910222),
    DistConstraint(9, 22, 1, 4.450789),
    DistConstraint(10, 22, 1, 4.145368),
    DistConstraint(11, 22, 1, 3.438209),
    DistConstraint(12, 22, 1, 3.906726),
    DistConstraint(13, 22, 1, 2.167323),
    DistConstraint(14, 22, 1, 2.851292),
    DistConstraint(15, 22, 1, 5.038529),
    DistConstraint(16, 22, 1, 4.986489),
    DistConstraint(17, 22, 1, 3.393839),
    DistConstraint(18, 22, 1, 1.080275),
    DistConstraint(19, 22, 1, 2.822729),
    DistConstraint(20, 22, 1, 4.266873),
    DistConstraint(21, 22, 1, 2.129777),
    DistConstraint(5, 23, 1, 6.437514),
    DistConstraint(9, 23, 1, 5.258844),
    DistConstraint(10, 23, 1, 5.571174),
    DistConstraint(11, 23, 1, 3.858649),
    DistConstraint(12, 23, 1, 3.386269),
    DistConstraint(13, 23, 1, 3.350602),
    DistConstraint(14, 23, 1, 4.607413),
    DistConstraint(15, 23, 1, 6.623072),
    DistConstraint(16, 23, 1, 4.269782),
    DistConstraint(17, 23, 1, 2.139060),
    DistConstraint(18, 23, 1, 2.125574),
    DistConstraint(19, 23, 1, 4.965133),
    DistConstraint(20, 23, 1, 2.441513),
    DistConstraint(21, 23, 1, 1.079013),
    DistConstraint(22, 23, 1, 2.456389),
    DistConstraint(1, 24, 1, 2.457634),
    DistConstraint(3, 24, 1, 1.524672),
    DistConstraint(4, 24, 1, 2.133947),
    DistConstraint(5, 24, 1, 2.501385),
    DistConstraint(26, 24, 1, 1.329000),
    DistConstraint(3, 25, 1, 2.400563),
    DistConstraint(24, 25, 1, 1.231053),
    DistConstraint(3, 26, 1, 2.424892),
    DistConstraint(24, 26, 1, 1.329000),
    DistConstraint(25, 26, 1, 2.250224),
    DistConstraint(28, 26, 1, 1.458081),
    DistConstraint(30, 26, 1, 2.452842),
    DistConstraint(40, 26, 1, 2.461671),
    DistConstraint(3, 27, 1, 2.545724),
    DistConstraint(24, 27, 1, 2.007267),
    DistConstraint(25, 27, 1, 3.120246),
    DistConstraint(26, 27, 1, 0.979986),
    DistConstraint(3, 28, 1, 3.803748),
    DistConstraint(24, 28, 1, 2.434942),
    DistConstraint(25, 28, 1, 2.773405),
    DistConstraint(26, 28, 1, 1.458081),
    DistConstraint(27, 28, 1, 2.107358),
    DistConstraint(40, 28, 1, 1.525000),
    DistConstraint(41, 28, 1, 2.400724),
    DistConstraint(26, 29, 1, 2.063952),
    DistConstraint(28, 29, 1, 1.080888),
    DistConstraint(26, 30, 1, 2.452842),
    DistConstraint(28, 30, 1, 1.531485),
    DistConstraint(29, 30, 1, 2.140267),
    DistConstraint(31, 30, 1, 1.080389),
    DistConstraint(33, 30, 1, 0.625257),
    DistConstraint(34, 30, 1, 1.516505),
    DistConstraint(40, 30, 1, 2.501819),
    DistConstraint(28, 31, 1, 2.144611),
    DistConstraint(30, 31, 1, 1.080389),
    DistConstraint(33, 31, 1, 0.880101),
    DistConstraint(34, 31, 1, 2.121281),
    DistConstraint(28, 32, 1, 2.143066),
    DistConstraint(30, 32, 1, 1.079571),
    DistConstraint(31, 32, 1, 1.761199),
    DistConstraint(28, 33, 1, 1.954231),
    DistConstraint(30, 33, 1, 0.625257),
    DistConstraint(31, 33, 1, 0.880101),
    DistConstraint(32, 33, 1, 0.881099),
    DistConstraint(28, 34, 1, 2.533282),
    DistConstraint(30, 34, 1, 1.516505),
    DistConstraint(31, 34, 1, 2.121281),
    DistConstraint(32, 34, 1, 2.121996),
    DistConstraint(33, 34, 1, 1.930722),
    DistConstraint(30, 35, 1, 2.394091),
    DistConstraint(34, 35, 1, 1.231098),
    DistConstraint(36, 35, 1, 2.245929),
    DistConstraint(37, 35, 1, 3.105081),
    DistConstraint(38, 35, 1, 2.461001),
    DistConstraint(30, 36, 1, 2.420160),
    DistConstraint(34, 36, 1, 1.327953),
    DistConstraint(35, 36, 1, 2.245929),
    DistConstraint(30, 37, 1, 2.515687),
    DistConstraint(34, 37, 1, 1.988646),
    DistConstraint(35, 37, 1, 3.105081),
    DistConstraint(36, 37, 1, 0.980184),
    DistConstraint(30, 38, 1, 3.319110),
    DistConstraint(34, 38, 1, 1.987932),
    DistConstraint(35, 38, 1, 2.461001),
    DistConstraint(36, 38, 1, 0.979957),
    DistConstraint(37, 38, 1, 1.727779),
    DistConstraint(30, 39, 1, 2.815366),
    DistConstraint(34, 39, 1, 1.790793),
    DistConstraint(35, 39, 1, 2.665069),
    DistConstraint(36, 39, 1, 0.462840),
    DistConstraint(37, 39, 1, 0.863918),
    DistConstraint(38, 39, 1, 0.863861),
    DistConstraint(26, 40, 1, 2.461671),
    DistConstraint(28, 40, 1, 1.525000),
    DistConstraint(29, 40, 1, 2.134832),
    DistConstraint(30, 40, 1, 2.501819),
    DistConstraint(28, 41, 1, 2.400724),
    DistConstraint(40, 41, 1, 1.231000),
};
static DistConstraints test_lo_bounds = {
    DistConstraint(2, 4, 1, 1.755000),
    DistConstraint(2, 5, 1, 2.295000),
    DistConstraint(1, 6, 1, 2.205000),
    DistConstraint(2, 6, 1, 1.755000),
    DistConstraint(4, 6, 1, 1.800000),
    DistConstraint(1, 7, 1, 2.205000),
    DistConstraint(2, 7, 1, 1.755000),
    DistConstraint(4, 7, 1, 1.800000),
    DistConstraint(1, 9, 1, 2.745000),
    DistConstraint(2, 9, 1, 2.295000),
    DistConstraint(4, 9, 1, 2.340000),
    DistConstraint(1, 10, 1, 2.745000),
    DistConstraint(2, 10, 1, 2.295000),
    DistConstraint(3, 10, 1, 2.880000),
    DistConstraint(4, 10, 1, 2.340000),
    DistConstraint(6, 10, 1, 2.340000),
    DistConstraint(7, 10, 1, 2.340000),
    DistConstraint(1, 11, 1, 2.745000),
    DistConstraint(2, 11, 1, 2.295000),
    DistConstraint(3, 11, 1, 2.880000),
    DistConstraint(4, 11, 1, 2.340000),
    DistConstraint(6, 11, 1, 2.340000),
    DistConstraint(7, 11, 1, 2.340000),
    DistConstraint(1, 12, 1, 2.745000),
    DistConstraint(2, 12, 1, 2.295000),
    DistConstraint(3, 12, 1, 2.880000),
    DistConstraint(4, 12, 1, 2.340000),
    DistConstraint(6, 12, 1, 2.340000),
    DistConstraint(7, 12, 1, 2.340000),
    DistConstraint(1, 13, 1, 2.745000),
    DistConstraint(2, 13, 1, 2.295000),
    DistConstraint(3, 13, 1, 2.880000),
    DistConstraint(4, 13, 1, 2.340000),
    DistConstraint(6, 13, 1, 2.340000),
    DistConstraint(7, 13, 1, 2.340000),
    DistConstraint(1, 14, 1, 2.610000),
    DistConstraint(2, 14, 1, 2.160000),
    DistConstraint(3, 14, 1, 2.745000),
    DistConstraint(4, 14, 1, 2.205000),
    DistConstraint(6, 14, 1, 2.205000),
    DistConstraint(7, 14, 1, 2.205000),
    DistConstraint(1, 15, 1, 2.205000),
    DistConstraint(2, 15, 1, 1.755000),
    DistConstraint(3, 15, 1, 2.340000),
    DistConstraint(4, 15, 1, 1.800000),
    DistConstraint(6, 15, 1, 1.800000),
    DistConstraint(7, 15, 1, 1.800000),
    DistConstraint(1, 16, 1, 2.205000),
    DistConstraint(2, 16, 1, 1.755000),
    DistConstraint(3, 16, 1, 2.340000),
    DistConstraint(4, 16, 1, 1.800000),
    DistConstraint(6, 16, 1, 1.800000),
    DistConstraint(7, 16, 1, 1.800000),
    DistConstraint(1, 17, 1, 2.745000),
    DistConstraint(2, 17, 1, 2.295000),
    DistConstraint(3, 17, 1, 2.880000),
    DistConstraint(4, 17, 1, 2.340000),
    DistConstraint(6, 17, 1, 2.340000),
    DistConstraint(7, 17, 1, 2.340000),
    DistConstraint(1, 18, 1, 2.745000),
    DistConstraint(2, 18, 1, 2.295000),
    DistConstraint(3, 18, 1, 2.880000),
    DistConstraint(4, 18, 1, 2.340000),
    DistConstraint(6, 18, 1, 2.340000),
    DistConstraint(7, 18, 1, 2.340000),
    DistConstraint(1, 19, 1, 2.160000),
    DistConstraint(2, 19, 1, 1.710000),
    DistConstraint(3, 19, 1, 2.295000),
    DistConstraint(4, 19, 1, 1.755000),
    DistConstraint(6, 19, 1, 1.755000),
    DistConstraint(7, 19, 1, 1.755000),
    DistConstraint(1, 20, 1, 2.205000),
    DistConstraint(2, 20, 1, 1.755000),
    DistConstraint(3, 20, 1, 2.340000),
    DistConstraint(4, 20, 1, 1.800000),
    DistConstraint(6, 20, 1, 1.800000),
    DistConstraint(7, 20, 1, 1.800000),
    DistConstraint(1, 21, 1, 2.745000),
    DistConstraint(2, 21, 1, 2.295000),
    DistConstraint(3, 21, 1, 2.880000),
    DistConstraint(4, 21, 1, 2.340000),
    DistConstraint(6, 21, 1, 2.340000),
    DistConstraint(7, 21, 1, 2.340000),
    DistConstraint(1, 22, 1, 2.205000),
    DistConstraint(2, 22, 1, 1.755000),
    DistConstraint(3, 22, 1, 2.340000),
    DistConstraint(4, 22, 1, 1.800000),
    DistConstraint(6, 22, 1, 1.800000),
    DistConstraint(7, 22, 1, 1.800000),
    DistConstraint(1, 23, 1, 2.205000),
    DistConstraint(2, 23, 1, 1.755000),
    DistConstraint(3, 23, 1, 2.340000),
    DistConstraint(4, 23, 1, 1.800000),
    DistConstraint(6, 23, 1, 1.800000),
    DistConstraint(7, 23, 1, 1.800000),
    DistConstraint(2, 24, 1, 2.205000),
    DistConstraint(6, 24, 1, 2.250000),
    DistConstraint(7, 24, 1, 2.250000),
    DistConstraint(9, 24, 1, 2.790000),
    DistConstraint(10, 24, 1, 2.790000),
    DistConstraint(11, 24, 1, 2.790000),
    DistConstraint(12, 24, 1, 2.790000),
    DistConstraint(13, 24, 1, 2.790000),
    DistConstraint(14, 24, 1, 2.655000),
    DistConstraint(15, 24, 1, 2.250000),
    DistConstraint(16, 24, 1, 2.250000),
    DistConstraint(17, 24, 1, 2.790000),
    DistConstraint(18, 24, 1, 2.790000),
    DistConstraint(19, 24, 1, 2.205000),
    DistConstraint(20, 24, 1, 2.250000),
    DistConstraint(21, 24, 1, 2.790000),
    DistConstraint(22, 24, 1, 2.250000),
    DistConstraint(23, 24, 1, 2.250000),
    DistConstraint(1, 25, 1, 2.475000),
    DistConstraint(2, 25, 1, 2.025000),
    DistConstraint(4, 25, 1, 2.070000),
    DistConstraint(5, 25, 1, 2.610000),
    DistConstraint(6, 25, 1, 2.070000),
    DistConstraint(7, 25, 1, 2.070000),
    DistConstraint(9, 25, 1, 2.610000),
    DistConstraint(10, 25, 1, 2.610000),
    DistConstraint(11, 25, 1, 2.610000),
    DistConstraint(12, 25, 1, 2.610000),
    DistConstraint(13, 25, 1, 2.610000),
    DistConstraint(14, 25, 1, 2.475000),
    DistConstraint(15, 25, 1, 2.070000),
    DistConstraint(16, 25, 1, 2.070000),
    DistConstraint(17, 25, 1, 2.610000),
    DistConstraint(18, 25, 1, 2.610000),
    DistConstraint(19, 25, 1, 2.025000),
    DistConstraint(20, 25, 1, 2.070000),
    DistConstraint(21, 25, 1, 2.610000),
    DistConstraint(22, 25, 1, 2.070000),
    DistConstraint(23, 25, 1, 2.070000),
    DistConstraint(1, 26, 1, 2.610000),
    DistConstraint(2, 26, 1, 2.160000),
    DistConstraint(4, 26, 1, 2.205000),
    DistConstraint(5, 26, 1, 2.745000),
    DistConstraint(6, 26, 1, 2.205000),
    DistConstraint(7, 26, 1, 2.205000),
    DistConstraint(9, 26, 1, 2.745000),
    DistConstraint(10, 26, 1, 2.745000),
    DistConstraint(11, 26, 1, 2.745000),
    DistConstraint(12, 26, 1, 2.745000),
    DistConstraint(13, 26, 1, 2.745000),
    DistConstraint(14, 26, 1, 2.610000),
    DistConstraint(15, 26, 1, 2.205000),
    DistConstraint(16, 26, 1, 2.205000),
    DistConstraint(17, 26, 1, 2.745000),
    DistConstraint(18, 26, 1, 2.745000),
    DistConstraint(19, 26, 1, 2.160000),
    DistConstraint(20, 26, 1, 2.205000),
    DistConstraint(21, 26, 1, 2.745000),
    DistConstraint(22, 26, 1, 2.205000),
    DistConstraint(23, 26, 1, 2.205000),
    DistConstraint(1, 27, 1, 2.160000),
    DistConstraint(2, 27, 1, 1.710000),
    DistConstraint(4, 27, 1, 1.755000),
    DistConstraint(5, 27, 1, 2.295000),
    DistConstraint(6, 27, 1, 1.755000),
    DistConstraint(7, 27, 1, 1.755000),
    DistConstraint(9, 27, 1, 2.295000),
    DistConstraint(10, 27, 1, 2.295000),
    DistConstraint(11, 27, 1, 2.295000),
    DistConstraint(12, 27, 1, 2.295000),
    DistConstraint(13, 27, 1, 2.295000),
    DistConstraint(14, 27, 1, 2.160000),
    DistConstraint(15, 27, 1, 1.755000),
    DistConstraint(16, 27, 1, 1.755000),
    DistConstraint(17, 27, 1, 2.295000),
    DistConstraint(18, 27, 1, 2.295000),
    DistConstraint(19, 27, 1, 1.710000),
    DistConstraint(20, 27, 1, 1.755000),
    DistConstraint(21, 27, 1, 2.295000),
    DistConstraint(22, 27, 1, 1.755000),
    DistConstraint(23, 27, 1, 1.755000),
    DistConstraint(1, 28, 1, 2.745000),
    DistConstraint(2, 28, 1, 2.295000),
    DistConstraint(4, 28, 1, 2.340000),
    DistConstraint(5, 28, 1, 2.880000),
    DistConstraint(6, 28, 1, 2.340000),
    DistConstraint(7, 28, 1, 2.340000),
    DistConstraint(9, 28, 1, 2.880000),
    DistConstraint(10, 28, 1, 2.880000),
    DistConstraint(11, 28, 1, 2.880000),
    DistConstraint(12, 28, 1, 2.880000),
    DistConstraint(13, 28, 1, 2.880000),
    DistConstraint(14, 28, 1, 2.745000),
    DistConstraint(15, 28, 1, 2.340000),
    DistConstraint(16, 28, 1, 2.340000),
    DistConstraint(17, 28, 1, 2.880000),
    DistConstraint(18, 28, 1, 2.880000),
    DistConstraint(19, 28, 1, 2.295000),
    DistConstraint(20, 28, 1, 2.340000),
    DistConstraint(21, 28, 1, 2.880000),
    DistConstraint(22, 28, 1, 2.340000),
    DistConstraint(23, 28, 1, 2.340000),
    DistConstraint(1, 29, 1, 2.205000),
    DistConstraint(2, 29, 1, 1.755000),
    DistConstraint(3, 29, 1, 2.340000),
    DistConstraint(4, 29, 1, 1.800000),
    DistConstraint(5, 29, 1, 2.340000),
    DistConstraint(6, 29, 1, 1.800000),
    DistConstraint(7, 29, 1, 1.800000),
    DistConstraint(9, 29, 1, 2.340000),
    DistConstraint(10, 29, 1, 2.340000),
    DistConstraint(11, 29, 1, 2.340000),
    DistConstraint(12, 29, 1, 2.340000),
    DistConstraint(13, 29, 1, 2.340000),
    DistConstraint(14, 29, 1, 2.205000),
    DistConstraint(15, 29, 1, 1.800000),
    DistConstraint(16, 29, 1, 1.800000),
    DistConstraint(17, 29, 1, 2.340000),
    DistConstraint(18, 29, 1, 2.340000),
    DistConstraint(19, 29, 1, 1.755000),
    DistConstraint(20, 29, 1, 1.800000),
    DistConstraint(21, 29, 1, 2.340000),
    DistConstraint(22, 29, 1, 1.800000),
    DistConstraint(23, 29, 1, 1.800000),
    DistConstraint(24, 29, 1, 2.250000),
    DistConstraint(25, 29, 1, 2.070000),
    DistConstraint(27, 29, 1, 1.755000),
    DistConstraint(1, 30, 1, 2.745000),
    DistConstraint(2, 30, 1, 2.295000),
    DistConstraint(3, 30, 1, 2.880000),
    DistConstraint(4, 30, 1, 2.340000),
    DistConstraint(5, 30, 1, 2.880000),
    DistConstraint(6, 30, 1, 2.340000),
    DistConstraint(7, 30, 1, 2.340000),
    DistConstraint(9, 30, 1, 2.880000),
    DistConstraint(10, 30, 1, 2.880000),
    DistConstraint(11, 30, 1, 2.880000),
    DistConstraint(12, 30, 1, 2.880000),
    DistConstraint(13, 30, 1, 2.880000),
    DistConstraint(14, 30, 1, 2.745000),
    DistConstraint(15, 30, 1, 2.340000),
    DistConstraint(16, 30, 1, 2.340000),
    DistConstraint(17, 30, 1, 2.880000),
    DistConstraint(18, 30, 1, 2.880000),
    DistConstraint(19, 30, 1, 2.295000),
    DistConstraint(20, 30, 1, 2.340000),
    DistConstraint(21, 30, 1, 2.880000),
    DistConstraint(22, 30, 1, 2.340000),
    DistConstraint(23, 30, 1, 2.340000),
    DistConstraint(24, 30, 1, 2.790000),
    DistConstraint(25, 30, 1, 2.610000),
    DistConstraint(27, 30, 1, 2.295000),
    DistConstraint(1, 31, 1, 2.205000),
    DistConstraint(2, 31, 1, 1.755000),
    DistConstraint(3, 31, 1, 2.340000),
    DistConstraint(4, 31, 1, 1.800000),
    DistConstraint(5, 31, 1, 2.340000),
    DistConstraint(6, 31, 1, 1.800000),
    DistConstraint(7, 31, 1, 1.800000),
    DistConstraint(9, 31, 1, 2.340000),
    DistConstraint(10, 31, 1, 2.340000),
    DistConstraint(11, 31, 1, 2.340000),
    DistConstraint(12, 31, 1, 2.340000),
    DistConstraint(13, 31, 1, 2.340000),
    DistConstraint(14, 31, 1, 2.205000),
    DistConstraint(15, 31, 1, 1.800000),
    DistConstraint(16, 31, 1, 1.800000),
    DistConstraint(17, 31, 1, 2.340000),
    DistConstraint(18, 31, 1, 2.340000),
    DistConstraint(19, 31, 1, 1.755000),
    DistConstraint(20, 31, 1, 1.800000),
    DistConstraint(21, 31, 1, 2.340000),
    DistConstraint(22, 31, 1, 1.800000),
    DistConstraint(23, 31, 1, 1.800000),
    DistConstraint(24, 31, 1, 2.250000),
    DistConstraint(25, 31, 1, 2.070000),
    DistConstraint(26, 31, 1, 2.205000),
    DistConstraint(27, 31, 1, 1.755000),
    DistConstraint(29, 31, 1, 1.800000),
    DistConstraint(1, 32, 1, 2.205000),
    DistConstraint(2, 32, 1, 1.755000),
    DistConstraint(3, 32, 1, 2.340000),
    DistConstraint(4, 32, 1, 1.800000),
    DistConstraint(5, 32, 1, 2.340000),
    DistConstraint(6, 32, 1, 1.800000),
    DistConstraint(7, 32, 1, 1.800000),
    DistConstraint(9, 32, 1, 2.340000),
    DistConstraint(10, 32, 1, 2.340000),
    DistConstraint(11, 32, 1, 2.340000),
    DistConstraint(12, 32, 1, 2.340000),
    DistConstraint(13, 32, 1, 2.340000),
    DistConstraint(14, 32, 1, 2.205000),
    DistConstraint(15, 32, 1, 1.800000),
    DistConstraint(16, 32, 1, 1.800000),
    DistConstraint(17, 32, 1, 2.340000),
    DistConstraint(18, 32, 1, 2.340000),
    DistConstraint(19, 32, 1, 1.755000),
    DistConstraint(20, 32, 1, 1.800000),
    DistConstraint(21, 32, 1, 2.340000),
    DistConstraint(22, 32, 1, 1.800000),
    DistConstraint(23, 32, 1, 1.800000),
    DistConstraint(24, 32, 1, 2.250000),
    DistConstraint(25, 32, 1, 2.070000),
    DistConstraint(26, 32, 1, 2.205000),
    DistConstraint(27, 32, 1, 1.755000),
    DistConstraint(29, 32, 1, 1.800000),
    DistConstraint(1, 34, 1, 2.655000),
    DistConstraint(2, 34, 1, 2.205000),
    DistConstraint(3, 34, 1, 2.790000),
    DistConstraint(4, 34, 1, 2.250000),
    DistConstraint(5, 34, 1, 2.790000),
    DistConstraint(6, 34, 1, 2.250000),
    DistConstraint(7, 34, 1, 2.250000),
    DistConstraint(9, 34, 1, 2.790000),
    DistConstraint(10, 34, 1, 2.790000),
    DistConstraint(11, 34, 1, 2.790000),
    DistConstraint(12, 34, 1, 2.790000),
    DistConstraint(13, 34, 1, 2.790000),
    DistConstraint(14, 34, 1, 2.655000),
    DistConstraint(15, 34, 1, 2.250000),
    DistConstraint(16, 34, 1, 2.250000),
    DistConstraint(17, 34, 1, 2.790000),
    DistConstraint(18, 34, 1, 2.790000),
    DistConstraint(19, 34, 1, 2.205000),
    DistConstraint(20, 34, 1, 2.250000),
    DistConstraint(21, 34, 1, 2.790000),
    DistConstraint(22, 34, 1, 2.250000),
    DistConstraint(23, 34, 1, 2.250000),
    DistConstraint(24, 34, 1, 2.700000),
    DistConstraint(25, 34, 1, 2.520000),
    DistConstraint(26, 34, 1, 2.655000),
    DistConstraint(27, 34, 1, 2.205000),
    DistConstraint(29, 34, 1, 2.250000),
    DistConstraint(1, 35, 1, 2.475000),
    DistConstraint(2, 35, 1, 2.025000),
    DistConstraint(3, 35, 1, 2.610000),
    DistConstraint(4, 35, 1, 2.070000),
    DistConstraint(5, 35, 1, 2.610000),
    DistConstraint(6, 35, 1, 2.070000),
    DistConstraint(7, 35, 1, 2.070000),
    DistConstraint(9, 35, 1, 2.610000),
    DistConstraint(10, 35, 1, 2.610000),
    DistConstraint(11, 35, 1, 2.610000),
    DistConstraint(12, 35, 1, 2.610000),
    DistConstraint(13, 35, 1, 2.610000),
    DistConstraint(14, 35, 1, 2.475000),
    DistConstraint(15, 35, 1, 2.070000),
    DistConstraint(16, 35, 1, 2.070000),
    DistConstraint(17, 35, 1, 2.610000),
    DistConstraint(18, 35, 1, 2.610000),
    DistConstraint(19, 35, 1, 2.025000),
    DistConstraint(20, 35, 1, 2.070000),
    DistConstraint(21, 35, 1, 2.610000),
    DistConstraint(22, 35, 1, 2.070000),
    DistConstraint(23, 35, 1, 2.070000),
    DistConstraint(24, 35, 1, 2.520000),
    DistConstraint(25, 35, 1, 2.340000),
    DistConstraint(26, 35, 1, 2.475000),
    DistConstraint(27, 35, 1, 2.025000),
    DistConstraint(28, 35, 1, 2.610000),
    DistConstraint(29, 35, 1, 2.070000),
    DistConstraint(31, 35, 1, 2.070000),
    DistConstraint(32, 35, 1, 2.070000),
    DistConstraint(1, 36, 1, 2.610000),
    DistConstraint(2, 36, 1, 2.160000),
    DistConstraint(3, 36, 1, 2.745000),
    DistConstraint(4, 36, 1, 2.205000),
    DistConstraint(5, 36, 1, 2.745000),
    DistConstraint(6, 36, 1, 2.205000),
    DistConstraint(7, 36, 1, 2.205000),
    DistConstraint(9, 36, 1, 2.745000),
    DistConstraint(10, 36, 1, 2.745000),
    DistConstraint(11, 36, 1, 2.745000),
    DistConstraint(12, 36, 1, 2.745000),
    DistConstraint(13, 36, 1, 2.745000),
    DistConstraint(14, 36, 1, 2.610000),
    DistConstraint(15, 36, 1, 2.205000),
    DistConstraint(16, 36, 1, 2.205000),
    DistConstraint(17, 36, 1, 2.745000),
    DistConstraint(18, 36, 1, 2.745000),
    DistConstraint(19, 36, 1, 2.160000),
    DistConstraint(20, 36, 1, 2.205000),
    DistConstraint(21, 36, 1, 2.745000),
    DistConstraint(22, 36, 1, 2.205000),
    DistConstraint(23, 36, 1, 2.205000),
    DistConstraint(24, 36, 1, 2.655000),
    DistConstraint(25, 36, 1, 2.475000),
    DistConstraint(26, 36, 1, 2.610000),
    DistConstraint(27, 36, 1, 2.160000),
    DistConstraint(28, 36, 1, 2.745000),
    DistConstraint(29, 36, 1, 2.205000),
    DistConstraint(31, 36, 1, 2.205000),
    DistConstraint(32, 36, 1, 2.205000),
    DistConstraint(1, 37, 1, 2.160000),
    DistConstraint(2, 37, 1, 1.710000),
    DistConstraint(3, 37, 1, 2.295000),
    DistConstraint(4, 37, 1, 1.755000),
    DistConstraint(5, 37, 1, 2.295000),
    DistConstraint(6, 37, 1, 1.755000),
    DistConstraint(7, 37, 1, 1.755000),
    DistConstraint(9, 37, 1, 2.295000),
    DistConstraint(10, 37, 1, 2.295000),
    DistConstraint(11, 37, 1, 2.295000),
    DistConstraint(12, 37, 1, 2.295000),
    DistConstraint(13, 37, 1, 2.295000),
    DistConstraint(14, 37, 1, 2.160000),
    DistConstraint(15, 37, 1, 1.755000),
    DistConstraint(16, 37, 1, 1.755000),
    DistConstraint(17, 37, 1, 2.295000),
    DistConstraint(18, 37, 1, 2.295000),
    DistConstraint(19, 37, 1, 1.710000),
    DistConstraint(20, 37, 1, 1.755000),
    DistConstraint(21, 37, 1, 2.295000),
    DistConstraint(22, 37, 1, 1.755000),
    DistConstraint(23, 37, 1, 1.755000),
    DistConstraint(24, 37, 1, 2.205000),
    DistConstraint(25, 37, 1, 2.025000),
    DistConstraint(26, 37, 1, 2.160000),
    DistConstraint(27, 37, 1, 1.710000),
    DistConstraint(28, 37, 1, 2.295000),
    DistConstraint(29, 37, 1, 1.755000),
    DistConstraint(31, 37, 1, 1.755000),
    DistConstraint(32, 37, 1, 1.755000),
    DistConstraint(1, 38, 1, 2.160000),
    DistConstraint(2, 38, 1, 1.710000),
    DistConstraint(3, 38, 1, 2.295000),
    DistConstraint(4, 38, 1, 1.755000),
    DistConstraint(5, 38, 1, 2.295000),
    DistConstraint(6, 38, 1, 1.755000),
    DistConstraint(7, 38, 1, 1.755000),
    DistConstraint(9, 38, 1, 2.295000),
    DistConstraint(10, 38, 1, 2.295000),
    DistConstraint(11, 38, 1, 2.295000),
    DistConstraint(12, 38, 1, 2.295000),
    DistConstraint(13, 38, 1, 2.295000),
    DistConstraint(14, 38, 1, 2.160000),
    DistConstraint(15, 38, 1, 1.755000),
    DistConstraint(16, 38, 1, 1.755000),
    DistConstraint(17, 38, 1, 2.295000),
    DistConstraint(18, 38, 1, 2.295000),
    DistConstraint(19, 38, 1, 1.710000),
    DistConstraint(20, 38, 1, 1.755000),
    DistConstraint(21, 38, 1, 2.295000),
    DistConstraint(22, 38, 1, 1.755000),
    DistConstraint(23, 38, 1, 1.755000),
    DistConstraint(24, 38, 1, 2.205000),
    DistConstraint(25, 38, 1, 2.025000),
    DistConstraint(26, 38, 1, 2.160000),
    DistConstraint(27, 38, 1, 1.710000),
    DistConstraint(28, 38, 1, 2.295000),
    DistConstraint(29, 38, 1, 1.755000),
    DistConstraint(31, 38, 1, 1.755000),
    DistConstraint(32, 38, 1, 1.755000),
    DistConstraint(1, 40, 1, 2.655000),
    DistConstraint(2, 40, 1, 2.205000),
    DistConstraint(3, 40, 1, 2.790000),
    DistConstraint(4, 40, 1, 2.250000),
    DistConstraint(5, 40, 1, 2.790000),
    DistConstraint(6, 40, 1, 2.250000),
    DistConstraint(7, 40, 1, 2.250000),
    DistConstraint(9, 40, 1, 2.790000),
    DistConstraint(10, 40, 1, 2.790000),
    DistConstraint(11, 40, 1, 2.790000),
    DistConstraint(12, 40, 1, 2.790000),
    DistConstraint(13, 40, 1, 2.790000),
    DistConstraint(14, 40, 1, 2.655000),
    DistConstraint(15, 40, 1, 2.250000),
    DistConstraint(16, 40, 1, 2.250000),
    DistConstraint(17, 40, 1, 2.790000),
    DistConstraint(18, 40, 1, 2.790000),
    DistConstraint(19, 40, 1, 2.205000),
    DistConstraint(20, 40, 1, 2.250000),
    DistConstraint(21, 40, 1, 2.790000),
    DistConstraint(22, 40, 1, 2.250000),
    DistConstraint(23, 40, 1, 2.250000),
    DistConstraint(24, 40, 1, 2.700000),
    DistConstraint(25, 40, 1, 2.520000),
    DistConstraint(27, 40, 1, 2.205000),
    DistConstraint(31, 40, 1, 2.250000),
    DistConstraint(32, 40, 1, 2.250000),
    DistConstraint(34, 40, 1, 2.700000),
    DistConstraint(35, 40, 1, 2.520000),
    DistConstraint(36, 40, 1, 2.655000),
    DistConstraint(37, 40, 1, 2.205000),
    DistConstraint(38, 40, 1, 2.205000),
    DistConstraint(1, 41, 1, 2.475000),
    DistConstraint(2, 41, 1, 2.025000),
    DistConstraint(3, 41, 1, 2.610000),
    DistConstraint(4, 41, 1, 2.070000),
    DistConstraint(5, 41, 1, 2.610000),
    DistConstraint(6, 41, 1, 2.070000),
    DistConstraint(7, 41, 1, 2.070000),
    DistConstraint(9, 41, 1, 2.610000),
    DistConstraint(10, 41, 1, 2.610000),
    DistConstraint(11, 41, 1, 2.610000),
    DistConstraint(12, 41, 1, 2.610000),
    DistConstraint(13, 41, 1, 2.610000),
    DistConstraint(14, 41, 1, 2.475000),
    DistConstraint(15, 41, 1, 2.070000),
    DistConstraint(16, 41, 1, 2.070000),
    DistConstraint(17, 41, 1, 2.610000),
    DistConstraint(18, 41, 1, 2.610000),
    DistConstraint(19, 41, 1, 2.025000),
    DistConstraint(20, 41, 1, 2.070000),
    DistConstraint(21, 41, 1, 2.610000),
    DistConstraint(22, 41, 1, 2.070000),
    DistConstraint(23, 41, 1, 2.070000),
    DistConstraint(24, 41, 1, 2.520000),
    DistConstraint(25, 41, 1, 2.340000),
    DistConstraint(26, 41, 1, 2.475000),
    DistConstraint(27, 41, 1, 2.025000),
    DistConstraint(29, 41, 1, 2.070000),
    DistConstraint(30, 41, 1, 2.610000),
    DistConstraint(31, 41, 1, 2.070000),
    DistConstraint(32, 41, 1, 2.070000),
    DistConstraint(34, 41, 1, 2.520000),
    DistConstraint(35, 41, 1, 2.340000),
    DistConstraint(36, 41, 1, 2.475000),
    DistConstraint(37, 41, 1, 2.025000),
    DistConstraint(38, 41, 1, 2.025000),
};
static DistConstraints test_up_bounds = {
    DistConstraint(4, 8, 0, 3.800000),
    DistConstraint(4, 2, 0, 3.300000),
    DistConstraint(8, 27, 0, 3.800000),
    DistConstraint(8, 2, 0, 3.200000),
    DistConstraint(2, 27, 0, 3.300000),
    DistConstraint(22, 23, 0, 2.700000),
    DistConstraint(20, 22, 0, 5.000000),
    DistConstraint(29, 33, 0, 2.700000),
    DistConstraint(29, 27, 0, 2.700000),
    DistConstraint(33, 27, 0, 2.700000),
    DistConstraint(27, 16, 0, 5.000000),
};
static double test_w[4] = { 2.000000, 1.000000, 1.000000, -1.000000 };
static double test_f[4] = { 10.000000, 10.000000, 10.000000, 10.000000 };
