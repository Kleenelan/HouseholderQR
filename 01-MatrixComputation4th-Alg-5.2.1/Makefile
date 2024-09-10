
#EXE := householder_vector
EXE += householder_qr_dec

all: $(EXE)

householder_vector: householder_vector.cpp utils.cpp
	g++ $^ -o $@ -lm -DBUILD_MAIN_VEC

householder_qr_dec: householder_qr_dec.cpp householder_vector.cpp utils.cpp
	g++ -g $^ -o $@ -lm -DBUILD_MAIN_DEC

.PHONY: clean
clean:
	-rm -rf $(EXE)
