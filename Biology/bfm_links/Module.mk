# svn $Id: Module.mk 429 2009-12-20 17:30:26Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
# Copyright (c) 2002-2010 The ROMS/TOMS Group             Kate Hedstrom :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

local_sub  := Waves/wwm_pgmcl_links

local_lib  := libWWM.a
local_src  := $(wildcard $(local_sub)/*.F)

$(eval $(call make-library,$(local_lib),$(local_src)))

$(eval $(compile-rules))

# swch = -unix -f95 -mpi
# @perl SWAN/switch.pl $(swch) SWAN/*.ftn


