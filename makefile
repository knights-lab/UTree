CC=gcc
.PHONY : all

all: utree-build utree-compress utree-search, utree-build_gg utree-search_gg

utree-build: itree.c ; $(CC) -m64 -std=gnu11 -Ofast -march=sandybridge itree.c -D BUILD PFBITS=24 IXTYPE=uint32_t

utree-compress: itree.c ; $(CC) -m64 -std=gnu11 -Ofast -march=sandybridge itree.c -D COMPRESS PFBITS=24 IXTYPE=uint32_t

utree-search: itree.c ; $(CC) -m64 -std=gnu11 -Ofast -march=sandybridge itree.c -D SEARCH PFBITS=24 IXTYPE=uint32_t

utree-build_gg: itree.c ; $(CC) -m64 -std=gnu11 -Ofast -march=sandybridge itree.c -D BUILD_GG PFBITS=24 IXTYPE=uint32_t

utree-search_gg: itree.c ; $(CC) -m64 -std=gnu11 -Ofast -march=sandybridge itree.c -D SEARCH_GG PFBITS=24 IXTYPE=uint32_t

