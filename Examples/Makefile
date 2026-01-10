include ../Makefile.common

test:
	./test_actions.exe
	./gauge.exe -hot -nc=3 -L=12x4x4x4 -output=gauge00.mdp
	./gauge.exe -input=gauge00.mdp -beta=6.0 -steps=10 -output=gauge01.mdp
	./meson.exe -input=gauge01.mdp -lkappa=0.1 -hkappa=0.1 -ps -clean
	./ks_meson.exe -input=gauge01.mdp -lmass=0.1 -hmass=0.1 -u0=1.0 -meson=5x5


