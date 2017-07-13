# Data_Nomenclature
Regular expressions for various sample type

## Reference of RNASeq sampleList 
#### "IRCR_LC14_342_T_RSq.1_fastqc.html", "IRCR_LC14_342_T_RSq.1_fastqc.zip", "IRCR_LC14_342_T_RSq.rpkm", "IRCR_LC14_346_RSq", "IRCR_LC15_446_RSq.html", "IRCR_BT16_1197_RSq", "IRCR_BT16_1290_T02_P0", "IRCR_BT16_1290_T02_P0_RSq"

### Regular expression
#### ^IRCR_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9]|_T0[1-9]_P0)_(RSq)+
<img src="https://regexper.com/#%5EIRCR_(%5BA-Z%5D%7B2%2C3%7D)(%5B0-9%5D%7B2%7D)_(%5B0-9%5D%7B3%2C4%7D)(%7C_T0%5B1-9%5D%7C_T0%5B1-9%5DP0)(RSq)%2B alt="alt text" title="Title" />

## Reference of Glioma Scan sampleList
#### "IRCR_GBM14_591_T01_BR1S1_P10_GS", "IRCR_BT16_1141_T_GS", "IRCR_BT16_1036_T04_GS", "IRCR_GBM14_593_P5_GS", "IRCR_GBM14_593_BR1_P0_GS", "HGF_IRCR_GBM13_210_BR1_P9_GS", "IRCR_GBM14_534_SP6_GS", "SNU_BT17_008_T_GS"
#### "S1311187_T_GS", "NS09_780_BR2_P25_GS", "S16124067_GS", "NS08_586_P33_GS", "NS08_586_BR1S7_P0_GS"

### Regular expression
#### ^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_B|_T|_T0[1-9])_(GS)+
#### ^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(GS)+
#### ^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(GS)+
<img src="https://regexper.com/#%5E(%7CHGF_)(IRCR%7CSNU%7CCHA%7CNCC%7CAMC)_(%5BA-Z%5D%7B2%2C3%7D)(%5B0-9%5D%7B2%7D)_(%5B0-9%5D%7B3%2C4%7D)(%7C_T0%5B1-9%5D)(_BR%5B1-9%5D%7C_BR%5B1-9%5DS%5B0-9%5D)(%7C_P%5B0-9%5D%7B1%2C2%7D%7C_SP%5B0-9%5D%7B1%2C2%7D)_(GS)%2B" alt="alt text" title="Title" />

#### ^S([0-9]{2})([0-9]{5,})(|_T)_(GS)+
#### ^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(GS)+
#### ^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(GS)+
