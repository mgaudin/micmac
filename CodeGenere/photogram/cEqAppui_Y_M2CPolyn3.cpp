// File Automatically generated by eLiSe
#include "general/all.h"
#include "private/all.h"
#include "cEqAppui_Y_M2CPolyn3.h"


cEqAppui_Y_M2CPolyn3::cEqAppui_Y_M2CPolyn3():
    cElCompiledFonc(1)
{
   AddIntRef (cIncIntervale("Intr",0,17));
   AddIntRef (cIncIntervale("Orient",17,23));
   Close(false);
}



void cEqAppui_Y_M2CPolyn3::ComputeVal()
{
   double tmp0_ = mCompCoord[17];
   double tmp1_ = mCompCoord[19];
   double tmp2_ = cos(tmp0_);
   double tmp3_ = cos(tmp1_);
   double tmp4_ = sin(tmp0_);
   double tmp5_ = mCompCoord[18];
   double tmp6_ = sin(tmp5_);
   double tmp7_ = -(tmp6_);
   double tmp8_ = sin(tmp1_);
   double tmp9_ = -(tmp4_);
   double tmp10_ = tmp2_*tmp7_;
   double tmp11_ = mCompCoord[20];
   double tmp12_ = mLocXTer-tmp11_;
   double tmp13_ = -(tmp8_);
   double tmp14_ = tmp4_*tmp7_;
   double tmp15_ = mCompCoord[21];
   double tmp16_ = mLocYTer-tmp15_;
   double tmp17_ = cos(tmp5_);
   double tmp18_ = mCompCoord[22];
   double tmp19_ = mLocZTer-tmp18_;
   double tmp20_ = mCompCoord[0];
   double tmp21_ = tmp9_*tmp13_;
   double tmp22_ = tmp10_*tmp3_;
   double tmp23_ = tmp21_+tmp22_;
   double tmp24_ = (tmp23_)*(tmp12_);
   double tmp25_ = tmp2_*tmp13_;
   double tmp26_ = tmp14_*tmp3_;
   double tmp27_ = tmp25_+tmp26_;
   double tmp28_ = (tmp27_)*(tmp16_);
   double tmp29_ = tmp24_+tmp28_;
   double tmp30_ = tmp17_*tmp3_;
   double tmp31_ = tmp30_*(tmp19_);
   double tmp32_ = tmp29_+tmp31_;
   double tmp33_ = mCompCoord[1];
   double tmp34_ = tmp2_*tmp17_;
   double tmp35_ = tmp34_*(tmp12_);
   double tmp36_ = tmp4_*tmp17_;
   double tmp37_ = tmp36_*(tmp16_);
   double tmp38_ = tmp35_+tmp37_;
   double tmp39_ = tmp6_*(tmp19_);
   double tmp40_ = tmp38_+tmp39_;
   double tmp41_ = (tmp40_)/(tmp32_);
   double tmp42_ = tmp20_*(tmp41_);
   double tmp43_ = tmp33_+tmp42_;
   double tmp44_ = (tmp43_)-mLocPolyn3_State_1_0;
   double tmp45_ = (tmp44_)/mLocPolyn3_State_0_0;
   double tmp46_ = mCompCoord[2];
   double tmp47_ = tmp9_*tmp3_;
   double tmp48_ = tmp10_*tmp8_;
   double tmp49_ = tmp47_+tmp48_;
   double tmp50_ = (tmp49_)*(tmp12_);
   double tmp51_ = tmp2_*tmp3_;
   double tmp52_ = tmp14_*tmp8_;
   double tmp53_ = tmp51_+tmp52_;
   double tmp54_ = (tmp53_)*(tmp16_);
   double tmp55_ = tmp50_+tmp54_;
   double tmp56_ = tmp17_*tmp8_;
   double tmp57_ = tmp56_*(tmp19_);
   double tmp58_ = tmp55_+tmp57_;
   double tmp59_ = (tmp58_)/(tmp32_);
   double tmp60_ = tmp20_*(tmp59_);
   double tmp61_ = tmp46_+tmp60_;
   double tmp62_ = (tmp61_)-mLocPolyn3_State_2_0;
   double tmp63_ = (tmp62_)/mLocPolyn3_State_0_0;
   double tmp64_ = (tmp45_)*(tmp45_);
   double tmp65_ = (tmp63_)*(tmp63_);

  mVal[0] = (mLocPolyn3_State_2_0+(((1-mCompCoord[3])*(tmp63_)+mCompCoord[4]*(tmp45_)+mCompCoord[5]*(tmp45_)*(tmp63_))-mCompCoord[6]*2*tmp65_+mCompCoord[8]*tmp64_)*mLocPolyn3_State_0_0+(mCompCoord[13]*tmp64_*(tmp45_)+mCompCoord[14]*(tmp63_)*(tmp45_)*(tmp45_)+mCompCoord[15]*tmp65_*(tmp45_)+mCompCoord[16]*(tmp63_)*tmp65_)*mLocPolyn3_State_0_0)-mLocYIm;

}


void cEqAppui_Y_M2CPolyn3::ComputeValDeriv()
{
   double tmp0_ = mCompCoord[17];
   double tmp1_ = mCompCoord[19];
   double tmp2_ = cos(tmp0_);
   double tmp3_ = cos(tmp1_);
   double tmp4_ = sin(tmp0_);
   double tmp5_ = mCompCoord[18];
   double tmp6_ = sin(tmp5_);
   double tmp7_ = -(tmp6_);
   double tmp8_ = sin(tmp1_);
   double tmp9_ = -(tmp4_);
   double tmp10_ = tmp2_*tmp7_;
   double tmp11_ = mCompCoord[20];
   double tmp12_ = mLocXTer-tmp11_;
   double tmp13_ = -(tmp8_);
   double tmp14_ = tmp4_*tmp7_;
   double tmp15_ = mCompCoord[21];
   double tmp16_ = mLocYTer-tmp15_;
   double tmp17_ = cos(tmp5_);
   double tmp18_ = mCompCoord[22];
   double tmp19_ = mLocZTer-tmp18_;
   double tmp20_ = mCompCoord[0];
   double tmp21_ = tmp9_*tmp13_;
   double tmp22_ = tmp10_*tmp3_;
   double tmp23_ = tmp21_+tmp22_;
   double tmp24_ = (tmp23_)*(tmp12_);
   double tmp25_ = tmp2_*tmp13_;
   double tmp26_ = tmp14_*tmp3_;
   double tmp27_ = tmp25_+tmp26_;
   double tmp28_ = (tmp27_)*(tmp16_);
   double tmp29_ = tmp24_+tmp28_;
   double tmp30_ = tmp17_*tmp3_;
   double tmp31_ = tmp30_*(tmp19_);
   double tmp32_ = tmp29_+tmp31_;
   double tmp33_ = mCompCoord[1];
   double tmp34_ = tmp2_*tmp17_;
   double tmp35_ = tmp34_*(tmp12_);
   double tmp36_ = tmp4_*tmp17_;
   double tmp37_ = tmp36_*(tmp16_);
   double tmp38_ = tmp35_+tmp37_;
   double tmp39_ = tmp6_*(tmp19_);
   double tmp40_ = tmp38_+tmp39_;
   double tmp41_ = (tmp40_)/(tmp32_);
   double tmp42_ = tmp20_*(tmp41_);
   double tmp43_ = tmp33_+tmp42_;
   double tmp44_ = (tmp43_)-mLocPolyn3_State_1_0;
   double tmp45_ = (tmp44_)/mLocPolyn3_State_0_0;
   double tmp46_ = mCompCoord[2];
   double tmp47_ = tmp9_*tmp3_;
   double tmp48_ = tmp10_*tmp8_;
   double tmp49_ = tmp47_+tmp48_;
   double tmp50_ = (tmp49_)*(tmp12_);
   double tmp51_ = tmp2_*tmp3_;
   double tmp52_ = tmp14_*tmp8_;
   double tmp53_ = tmp51_+tmp52_;
   double tmp54_ = (tmp53_)*(tmp16_);
   double tmp55_ = tmp50_+tmp54_;
   double tmp56_ = tmp17_*tmp8_;
   double tmp57_ = tmp56_*(tmp19_);
   double tmp58_ = tmp55_+tmp57_;
   double tmp59_ = (tmp58_)/(tmp32_);
   double tmp60_ = tmp20_*(tmp59_);
   double tmp61_ = tmp46_+tmp60_;
   double tmp62_ = (tmp61_)-mLocPolyn3_State_2_0;
   double tmp63_ = (tmp62_)/mLocPolyn3_State_0_0;
   double tmp64_ = (tmp45_)*(tmp45_);
   double tmp65_ = (tmp63_)*(tmp63_);
   double tmp66_ = mCompCoord[3];
   double tmp67_ = 1-tmp66_;
   double tmp68_ = ElSquare(mLocPolyn3_State_0_0);
   double tmp69_ = mCompCoord[4];
   double tmp70_ = (tmp41_)*mLocPolyn3_State_0_0;
   double tmp71_ = (tmp70_)/tmp68_;
   double tmp72_ = (tmp59_)*mLocPolyn3_State_0_0;
   double tmp73_ = (tmp72_)/tmp68_;
   double tmp74_ = mCompCoord[5];
   double tmp75_ = (tmp73_)*(tmp63_);
   double tmp76_ = mCompCoord[6];
   double tmp77_ = tmp76_*2;
   double tmp78_ = (tmp71_)*(tmp45_);
   double tmp79_ = mCompCoord[8];
   double tmp80_ = tmp78_+tmp78_;
   double tmp81_ = mCompCoord[13];
   double tmp82_ = (tmp73_)*(tmp45_);
   double tmp83_ = (tmp71_)*(tmp63_);
   double tmp84_ = (tmp63_)*(tmp45_);
   double tmp85_ = mCompCoord[14];
   double tmp86_ = tmp75_+tmp75_;
   double tmp87_ = mCompCoord[15];
   double tmp88_ = mCompCoord[16];
   double tmp89_ = mLocPolyn3_State_0_0/tmp68_;
   double tmp90_ = (tmp89_)*(tmp45_);
   double tmp91_ = tmp90_+tmp90_;
   double tmp92_ = (tmp89_)*(tmp63_);
   double tmp93_ = tmp92_+tmp92_;
   double tmp94_ = (tmp89_)*tmp65_;
   double tmp95_ = (tmp45_)*(tmp63_);
   double tmp96_ = tmp64_*(tmp45_);
   double tmp97_ = tmp84_*(tmp45_);
   double tmp98_ = tmp65_*(tmp45_);
   double tmp99_ = (tmp63_)*tmp65_;
   double tmp100_ = -(1);
   double tmp101_ = tmp100_*tmp4_;
   double tmp102_ = -(tmp2_);
   double tmp103_ = tmp101_*tmp7_;
   double tmp104_ = tmp102_*tmp13_;
   double tmp105_ = tmp103_*tmp3_;
   double tmp106_ = tmp104_+tmp105_;
   double tmp107_ = (tmp106_)*(tmp12_);
   double tmp108_ = tmp101_*tmp13_;
   double tmp109_ = tmp108_+tmp22_;
   double tmp110_ = (tmp109_)*(tmp16_);
   double tmp111_ = tmp107_+tmp110_;
   double tmp112_ = ElSquare(tmp32_);
   double tmp113_ = tmp101_*tmp17_;
   double tmp114_ = tmp113_*(tmp12_);
   double tmp115_ = tmp34_*(tmp16_);
   double tmp116_ = tmp114_+tmp115_;
   double tmp117_ = (tmp116_)*(tmp32_);
   double tmp118_ = (tmp40_)*(tmp111_);
   double tmp119_ = tmp117_-tmp118_;
   double tmp120_ = (tmp119_)/tmp112_;
   double tmp121_ = (tmp120_)*tmp20_;
   double tmp122_ = tmp121_*mLocPolyn3_State_0_0;
   double tmp123_ = (tmp122_)/tmp68_;
   double tmp124_ = tmp102_*tmp3_;
   double tmp125_ = tmp103_*tmp8_;
   double tmp126_ = tmp124_+tmp125_;
   double tmp127_ = (tmp126_)*(tmp12_);
   double tmp128_ = tmp101_*tmp3_;
   double tmp129_ = tmp128_+tmp48_;
   double tmp130_ = (tmp129_)*(tmp16_);
   double tmp131_ = tmp127_+tmp130_;
   double tmp132_ = (tmp131_)*(tmp32_);
   double tmp133_ = (tmp58_)*(tmp111_);
   double tmp134_ = tmp132_-tmp133_;
   double tmp135_ = (tmp134_)/tmp112_;
   double tmp136_ = (tmp135_)*tmp20_;
   double tmp137_ = tmp136_*mLocPolyn3_State_0_0;
   double tmp138_ = (tmp137_)/tmp68_;
   double tmp139_ = (tmp138_)*(tmp63_);
   double tmp140_ = (tmp123_)*(tmp45_);
   double tmp141_ = tmp140_+tmp140_;
   double tmp142_ = (tmp138_)*(tmp45_);
   double tmp143_ = (tmp123_)*(tmp63_);
   double tmp144_ = tmp139_+tmp139_;
   double tmp145_ = -(tmp17_);
   double tmp146_ = tmp145_*tmp2_;
   double tmp147_ = tmp145_*tmp4_;
   double tmp148_ = tmp100_*tmp6_;
   double tmp149_ = tmp146_*tmp3_;
   double tmp150_ = tmp149_*(tmp12_);
   double tmp151_ = tmp147_*tmp3_;
   double tmp152_ = tmp151_*(tmp16_);
   double tmp153_ = tmp150_+tmp152_;
   double tmp154_ = tmp148_*tmp3_;
   double tmp155_ = tmp154_*(tmp19_);
   double tmp156_ = tmp153_+tmp155_;
   double tmp157_ = tmp148_*tmp2_;
   double tmp158_ = tmp157_*(tmp12_);
   double tmp159_ = tmp148_*tmp4_;
   double tmp160_ = tmp159_*(tmp16_);
   double tmp161_ = tmp158_+tmp160_;
   double tmp162_ = tmp17_*(tmp19_);
   double tmp163_ = tmp161_+tmp162_;
   double tmp164_ = (tmp163_)*(tmp32_);
   double tmp165_ = (tmp40_)*(tmp156_);
   double tmp166_ = tmp164_-tmp165_;
   double tmp167_ = (tmp166_)/tmp112_;
   double tmp168_ = (tmp167_)*tmp20_;
   double tmp169_ = tmp168_*mLocPolyn3_State_0_0;
   double tmp170_ = (tmp169_)/tmp68_;
   double tmp171_ = tmp146_*tmp8_;
   double tmp172_ = tmp171_*(tmp12_);
   double tmp173_ = tmp147_*tmp8_;
   double tmp174_ = tmp173_*(tmp16_);
   double tmp175_ = tmp172_+tmp174_;
   double tmp176_ = tmp148_*tmp8_;
   double tmp177_ = tmp176_*(tmp19_);
   double tmp178_ = tmp175_+tmp177_;
   double tmp179_ = (tmp178_)*(tmp32_);
   double tmp180_ = (tmp58_)*(tmp156_);
   double tmp181_ = tmp179_-tmp180_;
   double tmp182_ = (tmp181_)/tmp112_;
   double tmp183_ = (tmp182_)*tmp20_;
   double tmp184_ = tmp183_*mLocPolyn3_State_0_0;
   double tmp185_ = (tmp184_)/tmp68_;
   double tmp186_ = (tmp185_)*(tmp63_);
   double tmp187_ = (tmp170_)*(tmp45_);
   double tmp188_ = tmp187_+tmp187_;
   double tmp189_ = (tmp185_)*(tmp45_);
   double tmp190_ = (tmp170_)*(tmp63_);
   double tmp191_ = tmp186_+tmp186_;
   double tmp192_ = tmp100_*tmp8_;
   double tmp193_ = -(tmp3_);
   double tmp194_ = tmp193_*tmp9_;
   double tmp195_ = tmp192_*tmp10_;
   double tmp196_ = tmp194_+tmp195_;
   double tmp197_ = (tmp196_)*(tmp12_);
   double tmp198_ = tmp193_*tmp2_;
   double tmp199_ = tmp192_*tmp14_;
   double tmp200_ = tmp198_+tmp199_;
   double tmp201_ = (tmp200_)*(tmp16_);
   double tmp202_ = tmp197_+tmp201_;
   double tmp203_ = tmp192_*tmp17_;
   double tmp204_ = tmp203_*(tmp19_);
   double tmp205_ = tmp202_+tmp204_;
   double tmp206_ = (tmp40_)*(tmp205_);
   double tmp207_ = -(tmp206_);
   double tmp208_ = tmp207_/tmp112_;
   double tmp209_ = (tmp208_)*tmp20_;
   double tmp210_ = tmp209_*mLocPolyn3_State_0_0;
   double tmp211_ = (tmp210_)/tmp68_;
   double tmp212_ = tmp192_*tmp9_;
   double tmp213_ = tmp3_*tmp10_;
   double tmp214_ = tmp212_+tmp213_;
   double tmp215_ = (tmp214_)*(tmp12_);
   double tmp216_ = tmp192_*tmp2_;
   double tmp217_ = tmp3_*tmp14_;
   double tmp218_ = tmp216_+tmp217_;
   double tmp219_ = (tmp218_)*(tmp16_);
   double tmp220_ = tmp215_+tmp219_;
   double tmp221_ = tmp3_*tmp17_;
   double tmp222_ = tmp221_*(tmp19_);
   double tmp223_ = tmp220_+tmp222_;
   double tmp224_ = (tmp223_)*(tmp32_);
   double tmp225_ = (tmp58_)*(tmp205_);
   double tmp226_ = tmp224_-tmp225_;
   double tmp227_ = (tmp226_)/tmp112_;
   double tmp228_ = (tmp227_)*tmp20_;
   double tmp229_ = tmp228_*mLocPolyn3_State_0_0;
   double tmp230_ = (tmp229_)/tmp68_;
   double tmp231_ = (tmp230_)*(tmp63_);
   double tmp232_ = (tmp211_)*(tmp45_);
   double tmp233_ = tmp232_+tmp232_;
   double tmp234_ = (tmp230_)*(tmp45_);
   double tmp235_ = (tmp211_)*(tmp63_);
   double tmp236_ = tmp231_+tmp231_;
   double tmp237_ = tmp100_*(tmp23_);
   double tmp238_ = tmp100_*tmp34_;
   double tmp239_ = tmp238_*(tmp32_);
   double tmp240_ = (tmp40_)*tmp237_;
   double tmp241_ = tmp239_-tmp240_;
   double tmp242_ = (tmp241_)/tmp112_;
   double tmp243_ = (tmp242_)*tmp20_;
   double tmp244_ = tmp243_*mLocPolyn3_State_0_0;
   double tmp245_ = (tmp244_)/tmp68_;
   double tmp246_ = tmp100_*(tmp49_);
   double tmp247_ = tmp246_*(tmp32_);
   double tmp248_ = (tmp58_)*tmp237_;
   double tmp249_ = tmp247_-tmp248_;
   double tmp250_ = (tmp249_)/tmp112_;
   double tmp251_ = (tmp250_)*tmp20_;
   double tmp252_ = tmp251_*mLocPolyn3_State_0_0;
   double tmp253_ = (tmp252_)/tmp68_;
   double tmp254_ = (tmp253_)*(tmp63_);
   double tmp255_ = (tmp245_)*(tmp45_);
   double tmp256_ = tmp255_+tmp255_;
   double tmp257_ = (tmp253_)*(tmp45_);
   double tmp258_ = (tmp245_)*(tmp63_);
   double tmp259_ = tmp254_+tmp254_;
   double tmp260_ = tmp100_*(tmp27_);
   double tmp261_ = tmp100_*tmp36_;
   double tmp262_ = tmp261_*(tmp32_);
   double tmp263_ = (tmp40_)*tmp260_;
   double tmp264_ = tmp262_-tmp263_;
   double tmp265_ = (tmp264_)/tmp112_;
   double tmp266_ = (tmp265_)*tmp20_;
   double tmp267_ = tmp266_*mLocPolyn3_State_0_0;
   double tmp268_ = (tmp267_)/tmp68_;
   double tmp269_ = tmp100_*(tmp53_);
   double tmp270_ = tmp269_*(tmp32_);
   double tmp271_ = (tmp58_)*tmp260_;
   double tmp272_ = tmp270_-tmp271_;
   double tmp273_ = (tmp272_)/tmp112_;
   double tmp274_ = (tmp273_)*tmp20_;
   double tmp275_ = tmp274_*mLocPolyn3_State_0_0;
   double tmp276_ = (tmp275_)/tmp68_;
   double tmp277_ = (tmp276_)*(tmp63_);
   double tmp278_ = (tmp268_)*(tmp45_);
   double tmp279_ = tmp278_+tmp278_;
   double tmp280_ = (tmp276_)*(tmp45_);
   double tmp281_ = (tmp268_)*(tmp63_);
   double tmp282_ = tmp277_+tmp277_;
   double tmp283_ = tmp100_*tmp30_;
   double tmp284_ = tmp148_*(tmp32_);
   double tmp285_ = (tmp40_)*tmp283_;
   double tmp286_ = tmp284_-tmp285_;
   double tmp287_ = (tmp286_)/tmp112_;
   double tmp288_ = (tmp287_)*tmp20_;
   double tmp289_ = tmp288_*mLocPolyn3_State_0_0;
   double tmp290_ = (tmp289_)/tmp68_;
   double tmp291_ = tmp100_*tmp56_;
   double tmp292_ = tmp291_*(tmp32_);
   double tmp293_ = (tmp58_)*tmp283_;
   double tmp294_ = tmp292_-tmp293_;
   double tmp295_ = (tmp294_)/tmp112_;
   double tmp296_ = (tmp295_)*tmp20_;
   double tmp297_ = tmp296_*mLocPolyn3_State_0_0;
   double tmp298_ = (tmp297_)/tmp68_;
   double tmp299_ = (tmp298_)*(tmp63_);
   double tmp300_ = (tmp290_)*(tmp45_);
   double tmp301_ = tmp300_+tmp300_;
   double tmp302_ = (tmp298_)*(tmp45_);
   double tmp303_ = (tmp290_)*(tmp63_);
   double tmp304_ = tmp299_+tmp299_;

  mVal[0] = (mLocPolyn3_State_2_0+(((tmp67_)*(tmp63_)+tmp69_*(tmp45_)+tmp74_*tmp95_)-tmp77_*tmp65_+tmp79_*tmp64_)*mLocPolyn3_State_0_0+(tmp81_*tmp96_+tmp85_*tmp97_+tmp87_*tmp98_+tmp88_*tmp99_)*mLocPolyn3_State_0_0)-mLocYIm;

  mCompDer[0][0] = (((tmp73_)*(tmp67_)+(tmp71_)*tmp69_+(tmp83_+tmp82_)*tmp74_)-(tmp86_)*tmp77_+(tmp80_)*tmp79_)*mLocPolyn3_State_0_0+(((tmp80_)*(tmp45_)+(tmp71_)*tmp64_)*tmp81_+((tmp82_+tmp83_)*(tmp45_)+(tmp71_)*tmp84_)*tmp85_+((tmp86_)*(tmp45_)+(tmp71_)*tmp65_)*tmp87_+((tmp73_)*tmp65_+(tmp86_)*(tmp63_))*tmp88_)*mLocPolyn3_State_0_0;
  mCompDer[0][1] = ((tmp89_)*tmp69_+tmp92_*tmp74_+(tmp91_)*tmp79_)*mLocPolyn3_State_0_0+(((tmp91_)*(tmp45_)+(tmp89_)*tmp64_)*tmp81_+(tmp92_*(tmp45_)+(tmp89_)*tmp84_)*tmp85_+tmp94_*tmp87_)*mLocPolyn3_State_0_0;
  mCompDer[0][2] = (((tmp89_)*(tmp67_)+tmp90_*tmp74_)-(tmp93_)*tmp77_)*mLocPolyn3_State_0_0+(tmp90_*(tmp45_)*tmp85_+(tmp93_)*(tmp45_)*tmp87_+(tmp94_+(tmp93_)*(tmp63_))*tmp88_)*mLocPolyn3_State_0_0;
  mCompDer[0][3] = tmp100_*(tmp63_)*mLocPolyn3_State_0_0;
  mCompDer[0][4] = (tmp45_)*mLocPolyn3_State_0_0;
  mCompDer[0][5] = tmp95_*mLocPolyn3_State_0_0;
  mCompDer[0][6] = -(2*tmp65_)*mLocPolyn3_State_0_0;
  mCompDer[0][7] = 0;
  mCompDer[0][8] = tmp64_*mLocPolyn3_State_0_0;
  mCompDer[0][9] = 0;
  mCompDer[0][10] = 0;
  mCompDer[0][11] = 0;
  mCompDer[0][12] = 0;
  mCompDer[0][13] = tmp96_*mLocPolyn3_State_0_0;
  mCompDer[0][14] = tmp97_*mLocPolyn3_State_0_0;
  mCompDer[0][15] = tmp98_*mLocPolyn3_State_0_0;
  mCompDer[0][16] = tmp99_*mLocPolyn3_State_0_0;
  mCompDer[0][17] = (((tmp138_)*(tmp67_)+(tmp123_)*tmp69_+(tmp143_+tmp142_)*tmp74_)-(tmp144_)*tmp77_+(tmp141_)*tmp79_)*mLocPolyn3_State_0_0+(((tmp141_)*(tmp45_)+(tmp123_)*tmp64_)*tmp81_+((tmp142_+tmp143_)*(tmp45_)+(tmp123_)*tmp84_)*tmp85_+((tmp144_)*(tmp45_)+(tmp123_)*tmp65_)*tmp87_+((tmp138_)*tmp65_+(tmp144_)*(tmp63_))*tmp88_)*mLocPolyn3_State_0_0;
  mCompDer[0][18] = (((tmp185_)*(tmp67_)+(tmp170_)*tmp69_+(tmp190_+tmp189_)*tmp74_)-(tmp191_)*tmp77_+(tmp188_)*tmp79_)*mLocPolyn3_State_0_0+(((tmp188_)*(tmp45_)+(tmp170_)*tmp64_)*tmp81_+((tmp189_+tmp190_)*(tmp45_)+(tmp170_)*tmp84_)*tmp85_+((tmp191_)*(tmp45_)+(tmp170_)*tmp65_)*tmp87_+((tmp185_)*tmp65_+(tmp191_)*(tmp63_))*tmp88_)*mLocPolyn3_State_0_0;
  mCompDer[0][19] = (((tmp230_)*(tmp67_)+(tmp211_)*tmp69_+(tmp235_+tmp234_)*tmp74_)-(tmp236_)*tmp77_+(tmp233_)*tmp79_)*mLocPolyn3_State_0_0+(((tmp233_)*(tmp45_)+(tmp211_)*tmp64_)*tmp81_+((tmp234_+tmp235_)*(tmp45_)+(tmp211_)*tmp84_)*tmp85_+((tmp236_)*(tmp45_)+(tmp211_)*tmp65_)*tmp87_+((tmp230_)*tmp65_+(tmp236_)*(tmp63_))*tmp88_)*mLocPolyn3_State_0_0;
  mCompDer[0][20] = (((tmp253_)*(tmp67_)+(tmp245_)*tmp69_+(tmp258_+tmp257_)*tmp74_)-(tmp259_)*tmp77_+(tmp256_)*tmp79_)*mLocPolyn3_State_0_0+(((tmp256_)*(tmp45_)+(tmp245_)*tmp64_)*tmp81_+((tmp257_+tmp258_)*(tmp45_)+(tmp245_)*tmp84_)*tmp85_+((tmp259_)*(tmp45_)+(tmp245_)*tmp65_)*tmp87_+((tmp253_)*tmp65_+(tmp259_)*(tmp63_))*tmp88_)*mLocPolyn3_State_0_0;
  mCompDer[0][21] = (((tmp276_)*(tmp67_)+(tmp268_)*tmp69_+(tmp281_+tmp280_)*tmp74_)-(tmp282_)*tmp77_+(tmp279_)*tmp79_)*mLocPolyn3_State_0_0+(((tmp279_)*(tmp45_)+(tmp268_)*tmp64_)*tmp81_+((tmp280_+tmp281_)*(tmp45_)+(tmp268_)*tmp84_)*tmp85_+((tmp282_)*(tmp45_)+(tmp268_)*tmp65_)*tmp87_+((tmp276_)*tmp65_+(tmp282_)*(tmp63_))*tmp88_)*mLocPolyn3_State_0_0;
  mCompDer[0][22] = (((tmp298_)*(tmp67_)+(tmp290_)*tmp69_+(tmp303_+tmp302_)*tmp74_)-(tmp304_)*tmp77_+(tmp301_)*tmp79_)*mLocPolyn3_State_0_0+(((tmp301_)*(tmp45_)+(tmp290_)*tmp64_)*tmp81_+((tmp302_+tmp303_)*(tmp45_)+(tmp290_)*tmp84_)*tmp85_+((tmp304_)*(tmp45_)+(tmp290_)*tmp65_)*tmp87_+((tmp298_)*tmp65_+(tmp304_)*(tmp63_))*tmp88_)*mLocPolyn3_State_0_0;
}


void cEqAppui_Y_M2CPolyn3::ComputeValDerivHessian()
{
  ELISE_ASSERT(false,"Foncteur cEqAppui_Y_M2CPolyn3 Has no Der Sec");
}

void cEqAppui_Y_M2CPolyn3::SetPolyn3_State_0_0(double aVal){ mLocPolyn3_State_0_0 = aVal;}
void cEqAppui_Y_M2CPolyn3::SetPolyn3_State_1_0(double aVal){ mLocPolyn3_State_1_0 = aVal;}
void cEqAppui_Y_M2CPolyn3::SetPolyn3_State_2_0(double aVal){ mLocPolyn3_State_2_0 = aVal;}
void cEqAppui_Y_M2CPolyn3::SetXTer(double aVal){ mLocXTer = aVal;}
void cEqAppui_Y_M2CPolyn3::SetYIm(double aVal){ mLocYIm = aVal;}
void cEqAppui_Y_M2CPolyn3::SetYTer(double aVal){ mLocYTer = aVal;}
void cEqAppui_Y_M2CPolyn3::SetZTer(double aVal){ mLocZTer = aVal;}



double * cEqAppui_Y_M2CPolyn3::AdrVarLocFromString(const std::string & aName)
{
   if (aName == "Polyn3_State_0_0") return & mLocPolyn3_State_0_0;
   if (aName == "Polyn3_State_1_0") return & mLocPolyn3_State_1_0;
   if (aName == "Polyn3_State_2_0") return & mLocPolyn3_State_2_0;
   if (aName == "XTer") return & mLocXTer;
   if (aName == "YIm") return & mLocYIm;
   if (aName == "YTer") return & mLocYTer;
   if (aName == "ZTer") return & mLocZTer;
   return 0;
}


cElCompiledFonc::cAutoAddEntry cEqAppui_Y_M2CPolyn3::mTheAuto("cEqAppui_Y_M2CPolyn3",cEqAppui_Y_M2CPolyn3::Alloc);


cElCompiledFonc *  cEqAppui_Y_M2CPolyn3::Alloc()
{  return new cEqAppui_Y_M2CPolyn3();
}


