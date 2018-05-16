/*
** svn $Id: hypoxia_srm_var.h 851 2017-06-23 23:22:54Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2017 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the Hypoxia Simple Respiration       **
**  Model biological variables that are used in input and output      **
**  NetCDF files. The metadata information is read from               **
**  "varinfo.dat".                                                    **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

/*
**  Model state biological tracers.
*/

            CASE ('idTvar(iOxyg)')
              idTvar(iOxyg)=varid
            CASE ('idTvar(iPO4_)')
              idTvar(iPO4_)=varid
            CASE ('idTvar(iNO3_)')
              idTvar(iNO3_)=varid
            CASE ('idTvar(iNH4_)')
              idTvar(iNH4_)=varid
            CASE ('idTvar(iO4n_)')
              idTvar(iO4n_)=varid
            CASE ('idTvar(iSiOH)')
              idTvar(iSiOH)=varid
            CASE ('idTvar(iN6r_)')
              idTvar(iN6r_)=varid

            CASE ('idTvar(iB1c_)')
              idTvar(iB1c_)=varid
            CASE ('idTvar(iB1n_)')
              idTvar(iB1n_)=varid
            CASE ('idTvar(iB1p_)')
              idTvar(iB1p_)=varid

            CASE ('idTvar(iP1c_)')
              idTvar(iP1c_)=varid
            CASE ('idTvar(iP1n_)')
              idTvar(iP1n_)=varid
            CASE ('idTvar(iP1p_)')
              idTvar(iP1p_)=varid
            CASE ('idTvar(iP1l_)')
              idTvar(iP1l_)=varid
            CASE ('idTvar(iP1s_)')
              idTvar(iP1s_)=varid

            CASE ('idTvar(iP2c_)')
              idTvar(iP2c_)=varid
            CASE ('idTvar(iP2n_)')
              idTvar(iP2n_)=varid
            CASE ('idTvar(iP2p_)')
              idTvar(iP2p_)=varid
            CASE ('idTvar(iP2l_)')
              idTvar(iP2l_)=varid

            CASE ('idTvar(iP3c_)')
              idTvar(iP3c_)=varid
            CASE ('idTvar(iP3n_)')
              idTvar(iP3n_)=varid
            CASE ('idTvar(iP3p_)')
              idTvar(iP3p_)=varid
            CASE ('idTvar(iP3l_)')
              idTvar(iP3l_)=varid

            CASE ('idTvar(iP4c_)')
              idTvar(iP4c_)=varid
            CASE ('idTvar(iP4n_)')
              idTvar(iP4n_)=varid
            CASE ('idTvar(iP4p_)')
              idTvar(iP4p_)=varid
            CASE ('idTvar(iP4l_)')
              idTvar(iP4l_)=varid

            CASE ('idTvar(iZ3c_)')
              idTvar(iZ3c_)=varid
            CASE ('idTvar(iZ3n_)')
              idTvar(iZ3n_)=varid
            CASE ('idTvar(iZ3p_)')
              idTvar(iZ3p_)=varid

            CASE ('idTvar(iZ4c_)')
              idTvar(iZ4c_)=varid
            CASE ('idTvar(iZ4n_)')
              idTvar(iZ4n_)=varid
            CASE ('idTvar(iZ4p_)')
              idTvar(iZ4p_)=varid

            CASE ('idTvar(iZ5c_)')
              idTvar(iZ5c_)=varid
            CASE ('idTvar(iZ5n_)')
              idTvar(iZ5n_)=varid
            CASE ('idTvar(iZ5p_)')
              idTvar(iZ5p_)=varid

            CASE ('idTvar(iZ6c_)')
              idTvar(iZ6c_)=varid
            CASE ('idTvar(iZ6n_)')
              idTvar(iZ6n_)=varid
            CASE ('idTvar(iZ6p_)')
              idTvar(iZ6p_)=varid

            CASE ('idTvar(iR1c_)')
              idTvar(iR1c_)=varid
            CASE ('idTvar(iR1n_)')
              idTvar(iR1n_)=varid
            CASE ('idTvar(iR1p_)')
              idTvar(iR1p_)=varid

            CASE ('idTvar(iR2c_)')
              idTvar(iR2c_)=varid

            CASE ('idTvar(iR3c_)')
              idTvar(iR3c_)=varid

            CASE ('idTvar(iR6c_)')
              idTvar(iR6c_)=varid
            CASE ('idTvar(iR6n_)')
              idTvar(iR6n_)=varid
            CASE ('idTvar(iR6p_)')
              idTvar(iR6p_)=varid
            CASE ('idTvar(iR6s_)')
              idTvar(iR6s_)=varid

            CASE ('idTvar(iO3c_)')
              idTvar(iO3c_)=varid
            CASE ('idTvar(iO3h_)')
              idTvar(iO3h_)=varid



            CASE ('idTvar(iEIR_)')
              idTvar(iEIR_)=varid
            CASE ('idTvar(iDIC_)')
              idTvar(iDIC_)=varid
            CASE ('idTvar(iChlo)')
              idTvar(iChlo)=varid

            CASE ('idTvar(siP1_)')
              idTvar(siP1_)=varid
            CASE ('idTvar(siP2_)')
              idTvar(siP2_)=varid
            CASE ('idTvar(siP3_)')
              idTvar(siP3_)=varid
            CASE ('idTvar(siP4_)')
              idTvar(siP4_)=varid
	      
            CASE ('idTvar(eiP1_)')
              idTvar(eiP1_)=varid
            CASE ('idTvar(eiP2_)')
              idTvar(eiP2_)=varid
            CASE ('idTvar(eiP3_)')
              idTvar(eiP3_)=varid
            CASE ('idTvar(eiP4_)')
              idTvar(eiP4_)=varid

            CASE ('idTvar(ruPTc)')
              idTvar(ruPTc)=varid
            CASE ('idTvar(ruZTc)')
              idTvar(ruZTc)=varid
            CASE ('idTvar(ixEPS)')
              idTvar(ixEPS)=varid
	      
/*
**  Biological tracers open boundary conditions.
*/

            CASE ('idTbry(iwest,iOxyg)')
              idTbry(iwest,iOxyg)=varid
            CASE ('idTbry(ieast,iOxyg)')
              idTbry(ieast,iOxyg)=varid
            CASE ('idTbry(isouth,iOxyg)')
              idTbry(isouth,iOxyg)=varid
            CASE ('idTbry(inorth,iOxyg)')
              idTbry(inorth,iOxyg)=varid

/*
**  Biological tracers point Source/Sinks (river runoff).
*/

            CASE ('idRtrc(iOxyg)')
              idRtrc(iOxyg)=varid
