CFLAGS = -g -std=c99 -I $(HOME)/include 
LDFLAGS = -g -L $(HOME)/lib -L /usr/X11R6/lib 
LIBS = -lm -lgraph -lX11

% : %.c

% : %.o $(OBJS)
	$(CC) $(LDFLAGS) $< $(LIBS) -o $@

%.o : %.c
	$(CC) -c $(CFLAGS)$<  -o $@

