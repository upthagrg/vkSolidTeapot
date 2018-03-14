sample.o:		sample.cpp
			g++ -std=gnu++11 -c -I.  sample.cpp


numbers.cpp:		sample.cpp
			rm -f numbers.cpp
			cat -n sample.cpp > numbers.cpp

sample-vert.spv:	sample-vert.vert
			glslangValidator -V sample-vert.vert  -o sample-vert.spv

sample-frag.spv:	sample-frag.frag
			glslangValidator -V sample-frag.frag  -o sample-frag.spv

shaders:		sample-vert.spv  sample-frag.spv


sample-vert-dis.txt:	sample-vert.vert
			rm -f sample-vert-dis.txt
			glslangValidator -H -V sample-vert.vert  > sample-vert-dis.txt

sample-frag-dis.txt:	sample-frag.frag
			rm -f sample-frag-dis.txt
			glslangValidator -H -V sample-frag.frag  > sample-frag-dis.txt

dis:			sample-vert-dis.txt  sample-frag-dis.txt


save:
			cp sample.cpp sample.save.cpp
			cp sample-vert.vert sample-vert.save.vert
			cp sample-frag.frag sample-frag.save.frag
