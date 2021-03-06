sampleList_RNASeq <- c("IRCR_LC14_342_T_RSq.1_fastqc.html", "IRCR_LC14_342_T_RSq.1_fastqc.zip", "IRCR_LC14_342_T_RSq.rpkm", "IRCR_LC14_346_RSq", "IRCR_LC15_446_RSq.html", "IRCR_BT16_1197_RSq", "IRCR_BT16_1290_T02_P0", "IRCR_BT16_1290_T02_P0_RSq")
grep("^IRCR_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9]|_T0[1-9]_P0)_(RSq)+", sampleList_RNASeq)


sampleList_GS_IRCR <- c("IRCR_GBM14_591_T01_BR1S1_P10_GS", "IRCR_BT16_1141_T_GS", "IRCR_BT16_1036_T04_GS", "IRCR_GBM14_593_P5_GS", "IRCR_GBM14_593_BR1_P0_GS", "HGF_IRCR_GBM13_210_BR1_P9_GS", "IRCR_GBM14_534_SP6_GS", "SNU_BT17_008_T_GS")

grep("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_B|_T|_T0[1-9])_(GS)+", sampleList_GS_IRCR)
grep("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(GS)+", sampleList_GS_IRCR)
grep("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(GS)+", sampleList_GS_IRCR)

sampleList_GS_OTHER <- c("S1311187_T_GS", "NS09_780_BR2_P25_GS", "S16124067_GS", "NS08_586_P33_GS", "NS08_586_BR1S7_P0_GS")
grep("^S([0-9]{2})([0-9]{5,})(|_T)_(GS)+", sampleList_GS_OTHER)
grep("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(GS)+", sampleList_GS_OTHER)
grep("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(GS)+", sampleList_GS_OTHER)
