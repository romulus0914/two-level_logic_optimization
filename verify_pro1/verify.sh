./verify2 case00
./verify2 case00.out
./verify2 case01
./verify2 case01.out
./verify2 case02
./verify2 case02.out
./verify2 case03
./verify2 case03.out
./verify2 case04
./verify2 case04.out
./verify2 case05
./verify2 case05.out
./verify2 case06
./verify2 case06.out
./verify2 case07
./verify2 case07.out

echo ====case00====
./verify1 case00.blif case00.out.blif
echo ====case01====
./verify1 case01.blif case01.out.blif
echo ====case02====
./verify1 case02.blif case02.out.blif
echo ====case03====
./verify1 case03.blif case03.out.blif
echo ====case04====
./verify1 case04.blif case04.out.blif
echo ====case05====
./verify1 case05.blif case05.out.blif
echo ====case06====
./verify1 case06.blif case06.out.blif
echo ====case07====
./verify1 case07.blif case07.out.blif

rm *.blif
