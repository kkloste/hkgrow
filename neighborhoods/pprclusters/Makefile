.PHONY : all clean rmake force_look

DIRS = hypercluster 

all: bin rmake
	-for d in $(DIRS); do (cp $$d/$$d bin/$$d ); done

bin: 
	mkdir -p bin
    
rmake: force_look
	-for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS) ); done
    
force_look :
	true
    
clean:
	$(RM) -rf bin
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

test: all
	./test_main.sh

