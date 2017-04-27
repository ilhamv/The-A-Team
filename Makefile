exec    = a.exe
cc      = g++
opt     = -O3
cflags  = -std=c++1y $(opt)
testdir = Unit_Testing

main    = Main.cpp
objects = $(patsubst %.cpp,%.o,$(filter-out $(main), $(wildcard *.cpp)))

.PHONY : all test clean

all :	$(objects) 
	@rm -f $(exec)
	@$(MAKE) $(exec) test

test :
#	@ cd $(testdir) && $(MAKE) # comment out this line to prevent unit testing (unit testing is slow)

%.o : %.cpp
	$(cc) $(cflags) -c $<

$(exec) : $(main)
	$(cc) $(cflags) $(objects) $< -o $@

clean :
	rm -f $(objects) $(exec)
	@ cd $(testdir) && $(MAKE) clean
