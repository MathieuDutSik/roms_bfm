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



            CASE ('EIR_')
              idTvar(iEIR_)=varid
            CASE ('iDIC_')
              idTvar(iDIC_)=varid
            CASE ('iChlo')
              idTvar(iChlo)=varid

            CASE ('siP1_')
              idTvar(siP1_)=varid
            CASE ('siP2_')
              idTvar(siP2_)=varid
            CASE ('siP3_')
              idTvar(siP3_)=varid
            CASE ('siP4_')
              idTvar(siP4_)=varid

            CASE ('eiP1_')
              idTvar(eiP1_)=varid
            CASE ('eiP2_')
              idTvar(eiP2_)=varid
            CASE ('eiP3_')
              idTvar(eiP3_)=varid
            CASE ('eiP4_')
              idTvar(eiP4_)=varid

            CASE ('ruPTc')
              idTvar(ruPTc)=varid
            CASE ('ruZTc')
              idTvar(ruZTc)=varid
            CASE ('ixEPS')
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
              Print *, 'Matching iOxyg, iOxyg=', iOxyg, ' varid=', varid
            CASE ('idRtrc(iPO4_)')
              idRtrc(iPO4_)=varid
              Print *, 'Matching iPO4 , iPO4_=', iPO4_, ' varid=', varid
            CASE ('idRtrc(iNO3_)')
              idRtrc(iNO3_)=varid
            CASE ('idRtrc(iNH4_)')
              idRtrc(iNH4_)=varid
            CASE ('idRtrc(iO4n_)')
              idRtrc(iO4n_)=varid
            CASE ('idRtrc(iSiOH)')
              idRtrc(iSiOH)=varid
            CASE ('idRtrc(iN6r_)')
              idRtrc(iN6r_)=varid

            CASE ('idRtrc(iB1c_)')
              idRtrc(iB1c_)=varid
            CASE ('idRtrc(iB1n_)')
              idRtrc(iB1n_)=varid
            CASE ('idRtrc(iB1p_)')
              idRtrc(iB1p_)=varid

            CASE ('idRtrc(iP1c_)')
              idRtrc(iP1c_)=varid
            CASE ('idRtrc(iP1n_)')
              idRtrc(iP1n_)=varid
            CASE ('idRtrc(iP1p_)')
              idRtrc(iP1p_)=varid
            CASE ('idRtrc(iP1l_)')
              idRtrc(iP1l_)=varid
            CASE ('idRtrc(iP1s_)')
              idRtrc(iP1s_)=varid

            CASE ('idRtrc(iP2c_)')
              idRtrc(iP2c_)=varid
            CASE ('idRtrc(iP2n_)')
              idRtrc(iP2n_)=varid
            CASE ('idRtrc(iP2p_)')
              idRtrc(iP2p_)=varid
            CASE ('idRtrc(iP2l_)')
              idRtrc(iP2l_)=varid

            CASE ('idRtrc(iP3c_)')
              idRtrc(iP3c_)=varid
            CASE ('idRtrc(iP3n_)')
              idRtrc(iP3n_)=varid
            CASE ('idRtrc(iP3p_)')
              idRtrc(iP3p_)=varid
            CASE ('idRtrc(iP3l_)')
              idRtrc(iP3l_)=varid

            CASE ('idRtrc(iP4c_)')
              idRtrc(iP4c_)=varid
            CASE ('idRtrc(iP4n_)')
              idRtrc(iP4n_)=varid
            CASE ('idRtrc(iP4p_)')
              idRtrc(iP4p_)=varid
            CASE ('idRtrc(iP4l_)')
              idRtrc(iP4l_)=varid

            CASE ('idRtrc(iZ3c_)')
              idRtrc(iZ3c_)=varid
            CASE ('idRtrc(iZ3n_)')
              idRtrc(iZ3n_)=varid
            CASE ('idRtrc(iZ3p_)')
              idRtrc(iZ3p_)=varid

            CASE ('idRtrc(iZ4c_)')
              idRtrc(iZ4c_)=varid
            CASE ('idRtrc(iZ4n_)')
              idRtrc(iZ4n_)=varid
            CASE ('idRtrc(iZ4p_)')
              idRtrc(iZ4p_)=varid

            CASE ('idRtrc(iZ5c_)')
              idRtrc(iZ5c_)=varid
            CASE ('idRtrc(iZ5n_)')
              idRtrc(iZ5n_)=varid
            CASE ('idRtrc(iZ5p_)')
              idRtrc(iZ5p_)=varid

            CASE ('idRtrc(iZ6c_)')
              idRtrc(iZ6c_)=varid
            CASE ('idRtrc(iZ6n_)')
              idRtrc(iZ6n_)=varid
            CASE ('idRtrc(iZ6p_)')
              idRtrc(iZ6p_)=varid

            CASE ('idRtrc(iR1c_)')
              idRtrc(iR1c_)=varid
            CASE ('idRtrc(iR1n_)')
              idRtrc(iR1n_)=varid
            CASE ('idRtrc(iR1p_)')
              idRtrc(iR1p_)=varid

            CASE ('idRtrc(iR2c_)')
              idRtrc(iR2c_)=varid

            CASE ('idRtrc(iR3c_)')
              idRtrc(iR3c_)=varid

            CASE ('idRtrc(iR6c_)')
              idRtrc(iR6c_)=varid
            CASE ('idRtrc(iR6n_)')
              idRtrc(iR6n_)=varid
            CASE ('idRtrc(iR6p_)')
              idRtrc(iR6p_)=varid
            CASE ('idRtrc(iR6s_)')
              idRtrc(iR6s_)=varid

            CASE ('idRtrc(iO3c_)')
              idRtrc(iO3c_)=varid
            CASE ('idRtrc(iO3h_)')
              idRtrc(iO3h_)=varid

            CASE ('idRtrc(iEIR_)')
              idRtrc(iEIR_)=varid
            CASE ('idRtrc(iDIC_)')
              idRtrc(iDIC_)=varid
            CASE ('idRtrc(iChlo)')
              idRtrc(iChlo)=varid

            CASE ('idRtrc(siP1_)')
              idRtrc(siP1_)=varid
            CASE ('idRtrc(siP2_)')
              idRtrc(siP2_)=varid
            CASE ('idRtrc(siP3_)')
              idRtrc(siP3_)=varid
            CASE ('idRtrc(siP4_)')
              idRtrc(siP4_)=varid

            CASE ('idRtrc(eiP1_)')
              idRtrc(eiP1_)=varid
            CASE ('idRtrc(eiP2_)')
              idRtrc(eiP2_)=varid
            CASE ('idRtrc(eiP3_)')
              idRtrc(eiP3_)=varid
            CASE ('idRtrc(eiP4_)')
              idRtrc(eiP4_)=varid

            CASE ('idRtrc(ruPTc)')
              idRtrc(ruPTc)=varid
            CASE ('idRtrc(ruZTc)')
              idRtrc(ruZTc)=varid
            CASE ('idRtrc(ixEPS)')
              idRtrc(ixEPS)=varid
