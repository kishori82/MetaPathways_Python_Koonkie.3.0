# The location of the expat directory
CC=gcc  
LEX=lex  
LEXFLAGS=-lfl
CFLAGS=-C

TARS= prodigal.v2_00.tar rpkm.tar trnascan.tar LAST-Plus.tar bwa-0.7.7.tar
SUBDIRS= prodigal.v2_00 rpkm trnascan  LAST-Plus bwa-0.7.7


build_folders:
	$(MAKE) $(CFLAGS) trnascan 
	$(MAKE) $(CFLAGS) rpkm 
	$(MAKE) $(CFLAGS) prodigal.v2_00 all
	$(MAKE) $(CFLAGS) LAST-Plus
	$(MAKE) $(CFLAGS) bwa-0.7.7

clean_folders:
	$(MAKE) $(CFLAGS) trnascan clean
	$(MAKE) $(CFLAGS) rpkm clean
	$(MAKE) $(CFLAGS) prodigal.v2_00 clean
	$(MAKE) $(CFLAGS) LAST-Plus clean
	$(MAKE) $(CFLAGS) bwa-0.7.7 clean

untar_folders: $(TARS)
	$(foreach tar, $(TARS), tar -xvf $(tar);)

tar_folders: $(SUBDIRS)
	$(foreach folder, $(SUBDIRS), tar -cvvf $(folder).tar $(folder);)

remove_folders:
	rm -rf $(SUBDIRS)

install_macosx:
	cp  trnascan/trnascan-1.4  ../macosx/
	cp  LAST-Plus/lastal+  ../macosx/ 
	cp  LAST-Plus/lastdb+  ../macosx/ 
	cp  bwa-0.7.7/bwa  ../macosx/ 
	cp  prodigal.v2_00/prodigal  ../macosx/ 
	cp  rpkm/rpkm  ../macosx/ 

remove_macosx:
	rm -rf  ../macosx/trnascan-1.4 
	rm -rf ../macosx/lastal+  
	rm -rf ../macosx/lastdb+  
	rm -rf ../macosx/bwa  
	rm -rf ../macosx/prodigal
	rm -rf ../macosx/rpkm 

install_redhat:
	cp  trnascan/trnascan-1.4  ../redhat/
	cp  LAST-Plus/lastal+  ../redhat/ 
	cp  LAST-Plus/lastdb+  ../redhat/ 
	cp  bwa-0.7.7/bwa  ../redhat/ 
	cp  prodigal.v2_00/prodigal  ../redhat/ 
	cp  rpkm/rpkm  ../redhat/ 

remove_redhat:
	rm -rf  ../redhat/trnascan-1.4 
	rm -rf ../redhat/lastal+  
	rm -rf ../redhat/lastdb+  
	rm -rf ../redhat/bwa  
	rm -rf ../redhat/prodigal
	rm -rf ../redhat/rpkm 

install_ubuntu:
	cp  trnascan/trnascan-1.4  ../ubuntu/
	cp  LAST-Plus/lastal+  ../ubuntu/ 
	cp  LAST-Plus/lastdb+  ../ubuntu/ 
	cp  bwa-0.7.7/bwa  ../ubuntu/ 
	cp  prodigal.v2_00/prodigal  ../ubuntu/ 
	cp  rpkm/rpkm  ../ubuntu/ 

remove_ubuntu:
	rm -rf  ../ubuntu/trnascan-1.4 
	rm -rf ../ubuntu/lastal+  
	rm -rf ../ubuntu/lastdb+  
	rm -rf ../ubuntu/bwa  
	rm -rf ../ubuntu/prodigal
	rm -rf ../ubuntu/rpkm 
