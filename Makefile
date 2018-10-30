psnw:ksw2_ggd.c cli.c ksw2.h
	$(CC) -Wall -g -O2 -o $@ ksw2_ggd.c cli.c

clean:
	rm -fr *.o *.dSYM psnw
