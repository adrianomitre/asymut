default : all

install : all
	cp asymut/asymut /usr/local/bin
	cp pda/pda /usr/local/bin
	cp peak_filter/peak_filter /usr/local/bin
	cp pitch2note/pitch2note /usr/local/bin
	cp keep_1st_ch.rb /usr/local/bin
	cp chrono_sort.rb /usr/local/bin

uninstall :
	rm /usr/local/bin/asymut
	rm /usr/local/bin/pda
	rm /usr/local/bin/peak_filter
	rm /usr/local/bin/pitch2note
	rm /usr/local/bin/keep_1st_ch.rb
	rm /usr/local/bin/chrono_sort.rb

all :
	make -C asymut
	make -C pda
	make -C peak_filter
	make -C pitch2note

clean :
	make clean -C asymut
	make clean -C pda
	make clean -C peak_filter
	make clean -C pitch2note
