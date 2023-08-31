RL:=$(shell lsb_release -si)

all:
ifeq ($(RL),"Ubuntu")
	make -f makefile.ubuntu
else ifeq ($(RL),"CentOS")
	make -f makefile.centos
else
	make -f makefile.osx
endif
