all: cover-figs neigh-figs nvsncp-figs core-figs pprgrow-figs codes data web

.phony: neigh-figs nvsncp-figs core-figs cover-figs pprgrow-figs data codes web

%.png : ../experiments/neighborhoods/%.eps
	convert -density 300 $< $@    
    
core-%.png : ../experiments/cores/core-%.eps
	convert -density 300 $< $@        
    
cover-%.png : ../experiments/overlaps/cover-%.eps
	convert -density 300 $< $@
	
pprgrow-%.png : ../experiments/grow/pprgrow-%.eps
	convert -density 300 $< $@	
    
graphs := itdk0304-cc  web-Google networkA-anonymized \
	rand-ff-25000-0.4 rand-ff-25000-0.49 \
	arxiv-ubc dblp-cc hollywood-2009-cc \
	ca-AstroPh-cc email-Enron-cc cond-mat-2005-fix-cc \
	Penn94 anony-interactions-oneyearA-cc \
	soc-LiveJournal1 oregon2_010526 \
	p2p-Gnutella25 as-22july06 itdk0304-cc \
	web-Google
	
neigh-comm.tar.gz : ../.git/index	
	git archive HEAD | gzip > neigh-comm.tar.gz

codes : neigh-comm.tar.gz

~/Dropbox/websites/purdue/codes/neighborhoods/neigh-comm-data.tar.gz : $(wildcard ../data/*.mat)
	tar czf ~/Dropbox/websites/purdue/codes/neighborhoods/neigh-comm-data.tar.gz $^ 

# $^ is the set of all prerequisites	
data : ~/Dropbox/websites/purdue/codes/neighborhoods/neigh-comm-data.tar.gz

index.md : ../README.md
	cp ../README.md index.md

web: index.md

deploy: all
	rsync -rtvp ./ ~/Dropbox/websites/purdue/codes/neighborhoods/ --exclude=Makefile --exclude=*~ --exclude=.DS_Store

neigh-figs: $(addprefix neigh-,$(addsuffix .png,$(graphs)))
nvsncp-figs: $(addprefix nvsncp-,$(addsuffix .png,$(graphs)))
core-figs: $(addprefix core-,$(addsuffix .png,$(graphs)))
cover-figs: $(addprefix cover-,$(addsuffix .png,$(graphs)))
pprgrow-figs: $(addprefix pprgrow-,$(addsuffix .png,$(graphs)))
