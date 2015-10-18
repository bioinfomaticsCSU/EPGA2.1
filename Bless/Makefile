PROG=bless
CC=g++
CFLAGS=-Wall
OPT=-O3
MAX_BOOST_ARRAY_SIZE=20
DEF=-D MAX_BOOST_ARRAY_SIZE=$(MAX_BOOST_ARRAY_SIZE)
SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) $(OBJS) -o $@

$(OBJS): %.o : %.cpp
	$(CC) $(CFLAGS) $(OPT) -c $(DEF) $<   

clean:
	rm -rf $(PROG) $(OBJS)
