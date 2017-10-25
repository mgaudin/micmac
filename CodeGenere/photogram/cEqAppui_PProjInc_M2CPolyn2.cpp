// File Automatically generated by eLiSe
#include "StdAfx.h"
#include "cEqAppui_PProjInc_M2CPolyn2.h"


cEqAppui_PProjInc_M2CPolyn2::cEqAppui_PProjInc_M2CPolyn2():
    cElCompiledFonc(2)
{
   AddIntRef (cIncIntervale("Intr",0,9));
   AddIntRef (cIncIntervale("Orient",9,15));
   AddIntRef (cIncIntervale("Tmp_PTer",15,18));
   Close(false);
}



void cEqAppui_PProjInc_M2CPolyn2::ComputeVal()
{
   double tmp0_ = mCompCoord[9];
   double tmp1_ = mCompCoord[10];
   double tmp2_ = cos(tmp1_);
   double tmp3_ = mCompCoord[15];
   double tmp4_ = mCompCoord[16];
   double tmp5_ = mCompCoord[17];
   double tmp6_ = sin(tmp0_);
   double tmp7_ = cos(tmp0_);
   double tmp8_ = sin(tmp1_);
   double tmp9_ = mCompCoord[11];
   double tmp10_ = mLocProjI_x * tmp3_;
   double tmp11_ = mLocProjP0_x + tmp10_;
   double tmp12_ = mLocProjJ_x * tmp4_;
   double tmp13_ = tmp11_ + tmp12_;
   double tmp14_ = mLocProjK_x * tmp5_;
   double tmp15_ = tmp13_ + tmp14_;
   double tmp16_ = mCompCoord[12];
   double tmp17_ = (tmp15_) - tmp16_;
   double tmp18_ = sin(tmp9_);
   double tmp19_ = -(tmp18_);
   double tmp20_ = -(tmp8_);
   double tmp21_ = cos(tmp9_);
   double tmp22_ = mLocProjI_y * tmp3_;
   double tmp23_ = mLocProjP0_y + tmp22_;
   double tmp24_ = mLocProjJ_y * tmp4_;
   double tmp25_ = tmp23_ + tmp24_;
   double tmp26_ = mLocProjK_y * tmp5_;
   double tmp27_ = tmp25_ + tmp26_;
   double tmp28_ = mCompCoord[13];
   double tmp29_ = (tmp27_) - tmp28_;
   double tmp30_ = mLocProjI_z * tmp3_;
   double tmp31_ = mLocProjP0_z + tmp30_;
   double tmp32_ = mLocProjJ_z * tmp4_;
   double tmp33_ = tmp31_ + tmp32_;
   double tmp34_ = mLocProjK_z * tmp5_;
   double tmp35_ = tmp33_ + tmp34_;
   double tmp36_ = mCompCoord[14];
   double tmp37_ = (tmp35_) - tmp36_;
   double tmp38_ = -(tmp6_);
   double tmp39_ = tmp7_ * tmp20_;
   double tmp40_ = tmp6_ * tmp20_;
   double tmp41_ = mCompCoord[0];
   double tmp42_ = tmp38_ * tmp19_;
   double tmp43_ = tmp39_ * tmp21_;
   double tmp44_ = tmp42_ + tmp43_;
   double tmp45_ = (tmp44_) * (tmp17_);
   double tmp46_ = tmp7_ * tmp19_;
   double tmp47_ = tmp40_ * tmp21_;
   double tmp48_ = tmp46_ + tmp47_;
   double tmp49_ = (tmp48_) * (tmp29_);
   double tmp50_ = tmp45_ + tmp49_;
   double tmp51_ = tmp2_ * tmp21_;
   double tmp52_ = tmp51_ * (tmp37_);
   double tmp53_ = tmp50_ + tmp52_;
   double tmp54_ = tmp41_ / (tmp53_);
   double tmp55_ = tmp7_ * tmp2_;
   double tmp56_ = tmp55_ * (tmp17_);
   double tmp57_ = tmp6_ * tmp2_;
   double tmp58_ = tmp57_ * (tmp29_);
   double tmp59_ = tmp56_ + tmp58_;
   double tmp60_ = tmp8_ * (tmp37_);
   double tmp61_ = tmp59_ + tmp60_;
   double tmp62_ = (tmp61_) * (tmp54_);
   double tmp63_ = mCompCoord[1];
   double tmp64_ = tmp62_ + tmp63_;
   double tmp65_ = (tmp64_) - mLocPolyn2_State_1_0;
   double tmp66_ = (tmp65_) / mLocPolyn2_State_0_0;
   double tmp67_ = tmp38_ * tmp21_;
   double tmp68_ = tmp39_ * tmp18_;
   double tmp69_ = tmp67_ + tmp68_;
   double tmp70_ = (tmp69_) * (tmp17_);
   double tmp71_ = tmp7_ * tmp21_;
   double tmp72_ = tmp40_ * tmp18_;
   double tmp73_ = tmp71_ + tmp72_;
   double tmp74_ = (tmp73_) * (tmp29_);
   double tmp75_ = tmp70_ + tmp74_;
   double tmp76_ = tmp2_ * tmp18_;
   double tmp77_ = tmp76_ * (tmp37_);
   double tmp78_ = tmp75_ + tmp77_;
   double tmp79_ = (tmp78_) * (tmp54_);
   double tmp80_ = mCompCoord[2];
   double tmp81_ = tmp79_ + tmp80_;
   double tmp82_ = (tmp81_) - mLocPolyn2_State_2_0;
   double tmp83_ = (tmp82_) / mLocPolyn2_State_0_0;
   double tmp84_ = mCompCoord[3];
   double tmp85_ = mCompCoord[4];
   double tmp86_ = mCompCoord[5];
   double tmp87_ = (tmp66_) * (tmp83_);
   double tmp88_ = mCompCoord[6];
   double tmp89_ = (tmp83_) * (tmp83_);
   double tmp90_ = (tmp66_) * (tmp66_);

  mVal[0] = ((mLocPolyn2_State_1_0 + (((1 + tmp84_) * (tmp66_) + tmp85_ * (tmp83_)) - tmp86_ * 2 * tmp90_ + tmp88_ * tmp87_ + mCompCoord[7] * tmp89_) * mLocPolyn2_State_0_0) - mLocXIm) * mLocScNorm;

  mVal[1] = ((mLocPolyn2_State_2_0 + (((1 - tmp84_) * (tmp83_) + tmp85_ * (tmp66_) + tmp86_ * tmp87_) - tmp88_ * 2 * tmp89_ + mCompCoord[8] * tmp90_) * mLocPolyn2_State_0_0) - mLocYIm) * mLocScNorm;

}


void cEqAppui_PProjInc_M2CPolyn2::ComputeValDeriv()
{
   double tmp0_ = mCompCoord[9];
   double tmp1_ = mCompCoord[10];
   double tmp2_ = cos(tmp1_);
   double tmp3_ = mCompCoord[15];
   double tmp4_ = mCompCoord[16];
   double tmp5_ = mCompCoord[17];
   double tmp6_ = sin(tmp0_);
   double tmp7_ = cos(tmp0_);
   double tmp8_ = sin(tmp1_);
   double tmp9_ = mCompCoord[11];
   double tmp10_ = mLocProjI_x * tmp3_;
   double tmp11_ = mLocProjP0_x + tmp10_;
   double tmp12_ = mLocProjJ_x * tmp4_;
   double tmp13_ = tmp11_ + tmp12_;
   double tmp14_ = mLocProjK_x * tmp5_;
   double tmp15_ = tmp13_ + tmp14_;
   double tmp16_ = mCompCoord[12];
   double tmp17_ = (tmp15_) - tmp16_;
   double tmp18_ = sin(tmp9_);
   double tmp19_ = -(tmp18_);
   double tmp20_ = -(tmp8_);
   double tmp21_ = cos(tmp9_);
   double tmp22_ = mLocProjI_y * tmp3_;
   double tmp23_ = mLocProjP0_y + tmp22_;
   double tmp24_ = mLocProjJ_y * tmp4_;
   double tmp25_ = tmp23_ + tmp24_;
   double tmp26_ = mLocProjK_y * tmp5_;
   double tmp27_ = tmp25_ + tmp26_;
   double tmp28_ = mCompCoord[13];
   double tmp29_ = (tmp27_) - tmp28_;
   double tmp30_ = mLocProjI_z * tmp3_;
   double tmp31_ = mLocProjP0_z + tmp30_;
   double tmp32_ = mLocProjJ_z * tmp4_;
   double tmp33_ = tmp31_ + tmp32_;
   double tmp34_ = mLocProjK_z * tmp5_;
   double tmp35_ = tmp33_ + tmp34_;
   double tmp36_ = mCompCoord[14];
   double tmp37_ = (tmp35_) - tmp36_;
   double tmp38_ = -(tmp6_);
   double tmp39_ = tmp7_ * tmp20_;
   double tmp40_ = tmp6_ * tmp20_;
   double tmp41_ = mCompCoord[0];
   double tmp42_ = tmp38_ * tmp19_;
   double tmp43_ = tmp39_ * tmp21_;
   double tmp44_ = tmp42_ + tmp43_;
   double tmp45_ = (tmp44_) * (tmp17_);
   double tmp46_ = tmp7_ * tmp19_;
   double tmp47_ = tmp40_ * tmp21_;
   double tmp48_ = tmp46_ + tmp47_;
   double tmp49_ = (tmp48_) * (tmp29_);
   double tmp50_ = tmp45_ + tmp49_;
   double tmp51_ = tmp2_ * tmp21_;
   double tmp52_ = tmp51_ * (tmp37_);
   double tmp53_ = tmp50_ + tmp52_;
   double tmp54_ = tmp41_ / (tmp53_);
   double tmp55_ = tmp7_ * tmp2_;
   double tmp56_ = tmp55_ * (tmp17_);
   double tmp57_ = tmp6_ * tmp2_;
   double tmp58_ = tmp57_ * (tmp29_);
   double tmp59_ = tmp56_ + tmp58_;
   double tmp60_ = tmp8_ * (tmp37_);
   double tmp61_ = tmp59_ + tmp60_;
   double tmp62_ = (tmp61_) * (tmp54_);
   double tmp63_ = mCompCoord[1];
   double tmp64_ = tmp62_ + tmp63_;
   double tmp65_ = (tmp64_) - mLocPolyn2_State_1_0;
   double tmp66_ = (tmp65_) / mLocPolyn2_State_0_0;
   double tmp67_ = tmp38_ * tmp21_;
   double tmp68_ = tmp39_ * tmp18_;
   double tmp69_ = tmp67_ + tmp68_;
   double tmp70_ = (tmp69_) * (tmp17_);
   double tmp71_ = tmp7_ * tmp21_;
   double tmp72_ = tmp40_ * tmp18_;
   double tmp73_ = tmp71_ + tmp72_;
   double tmp74_ = (tmp73_) * (tmp29_);
   double tmp75_ = tmp70_ + tmp74_;
   double tmp76_ = tmp2_ * tmp18_;
   double tmp77_ = tmp76_ * (tmp37_);
   double tmp78_ = tmp75_ + tmp77_;
   double tmp79_ = (tmp78_) * (tmp54_);
   double tmp80_ = mCompCoord[2];
   double tmp81_ = tmp79_ + tmp80_;
   double tmp82_ = (tmp81_) - mLocPolyn2_State_2_0;
   double tmp83_ = (tmp82_) / mLocPolyn2_State_0_0;
   double tmp84_ = mCompCoord[3];
   double tmp85_ = 1 + tmp84_;
   double tmp86_ = ElSquare(tmp53_);
   double tmp87_ = (tmp53_) / tmp86_;
   double tmp88_ = ElSquare(mLocPolyn2_State_0_0);
   double tmp89_ = mCompCoord[4];
   double tmp90_ = (tmp87_) * (tmp61_);
   double tmp91_ = tmp90_ * mLocPolyn2_State_0_0;
   double tmp92_ = (tmp91_) / tmp88_;
   double tmp93_ = (tmp92_) * (tmp66_);
   double tmp94_ = mCompCoord[5];
   double tmp95_ = tmp94_ * 2;
   double tmp96_ = (tmp87_) * (tmp78_);
   double tmp97_ = tmp96_ * mLocPolyn2_State_0_0;
   double tmp98_ = (tmp97_) / tmp88_;
   double tmp99_ = mCompCoord[6];
   double tmp100_ = (tmp98_) * (tmp83_);
   double tmp101_ = mCompCoord[7];
   double tmp102_ = mLocPolyn2_State_0_0 / tmp88_;
   double tmp103_ = (tmp102_) * (tmp66_);
   double tmp104_ = (tmp102_) * (tmp83_);
   double tmp105_ = (tmp66_) * (tmp66_);
   double tmp106_ = (tmp66_) * (tmp83_);
   double tmp107_ = (tmp83_) * (tmp83_);
   double tmp108_ = -(1);
   double tmp109_ = tmp108_ * tmp6_;
   double tmp110_ = -(tmp7_);
   double tmp111_ = tmp109_ * tmp20_;
   double tmp112_ = tmp110_ * tmp19_;
   double tmp113_ = tmp111_ * tmp21_;
   double tmp114_ = tmp112_ + tmp113_;
   double tmp115_ = (tmp114_) * (tmp17_);
   double tmp116_ = tmp109_ * tmp19_;
   double tmp117_ = tmp116_ + tmp43_;
   double tmp118_ = (tmp117_) * (tmp29_);
   double tmp119_ = tmp115_ + tmp118_;
   double tmp120_ = tmp41_ * (tmp119_);
   double tmp121_ = -(tmp120_);
   double tmp122_ = tmp121_ / tmp86_;
   double tmp123_ = tmp109_ * tmp2_;
   double tmp124_ = tmp123_ * (tmp17_);
   double tmp125_ = tmp55_ * (tmp29_);
   double tmp126_ = tmp124_ + tmp125_;
   double tmp127_ = (tmp126_) * (tmp54_);
   double tmp128_ = (tmp122_) * (tmp61_);
   double tmp129_ = tmp127_ + tmp128_;
   double tmp130_ = (tmp129_) * mLocPolyn2_State_0_0;
   double tmp131_ = (tmp130_) / tmp88_;
   double tmp132_ = (tmp131_) * (tmp66_);
   double tmp133_ = tmp110_ * tmp21_;
   double tmp134_ = tmp111_ * tmp18_;
   double tmp135_ = tmp133_ + tmp134_;
   double tmp136_ = (tmp135_) * (tmp17_);
   double tmp137_ = tmp109_ * tmp21_;
   double tmp138_ = tmp137_ + tmp68_;
   double tmp139_ = (tmp138_) * (tmp29_);
   double tmp140_ = tmp136_ + tmp139_;
   double tmp141_ = (tmp140_) * (tmp54_);
   double tmp142_ = (tmp122_) * (tmp78_);
   double tmp143_ = tmp141_ + tmp142_;
   double tmp144_ = (tmp143_) * mLocPolyn2_State_0_0;
   double tmp145_ = (tmp144_) / tmp88_;
   double tmp146_ = (tmp145_) * (tmp83_);
   double tmp147_ = tmp108_ * tmp8_;
   double tmp148_ = -(tmp2_);
   double tmp149_ = tmp148_ * tmp7_;
   double tmp150_ = tmp148_ * tmp6_;
   double tmp151_ = tmp149_ * tmp21_;
   double tmp152_ = tmp151_ * (tmp17_);
   double tmp153_ = tmp150_ * tmp21_;
   double tmp154_ = tmp153_ * (tmp29_);
   double tmp155_ = tmp152_ + tmp154_;
   double tmp156_ = tmp147_ * tmp21_;
   double tmp157_ = tmp156_ * (tmp37_);
   double tmp158_ = tmp155_ + tmp157_;
   double tmp159_ = tmp41_ * (tmp158_);
   double tmp160_ = -(tmp159_);
   double tmp161_ = tmp160_ / tmp86_;
   double tmp162_ = tmp147_ * tmp7_;
   double tmp163_ = tmp162_ * (tmp17_);
   double tmp164_ = tmp147_ * tmp6_;
   double tmp165_ = tmp164_ * (tmp29_);
   double tmp166_ = tmp163_ + tmp165_;
   double tmp167_ = tmp2_ * (tmp37_);
   double tmp168_ = tmp166_ + tmp167_;
   double tmp169_ = (tmp168_) * (tmp54_);
   double tmp170_ = (tmp161_) * (tmp61_);
   double tmp171_ = tmp169_ + tmp170_;
   double tmp172_ = (tmp171_) * mLocPolyn2_State_0_0;
   double tmp173_ = (tmp172_) / tmp88_;
   double tmp174_ = (tmp173_) * (tmp66_);
   double tmp175_ = tmp149_ * tmp18_;
   double tmp176_ = tmp175_ * (tmp17_);
   double tmp177_ = tmp150_ * tmp18_;
   double tmp178_ = tmp177_ * (tmp29_);
   double tmp179_ = tmp176_ + tmp178_;
   double tmp180_ = tmp147_ * tmp18_;
   double tmp181_ = tmp180_ * (tmp37_);
   double tmp182_ = tmp179_ + tmp181_;
   double tmp183_ = (tmp182_) * (tmp54_);
   double tmp184_ = (tmp161_) * (tmp78_);
   double tmp185_ = tmp183_ + tmp184_;
   double tmp186_ = (tmp185_) * mLocPolyn2_State_0_0;
   double tmp187_ = (tmp186_) / tmp88_;
   double tmp188_ = (tmp187_) * (tmp83_);
   double tmp189_ = -(tmp21_);
   double tmp190_ = tmp108_ * tmp18_;
   double tmp191_ = tmp189_ * tmp38_;
   double tmp192_ = tmp190_ * tmp39_;
   double tmp193_ = tmp191_ + tmp192_;
   double tmp194_ = (tmp193_) * (tmp17_);
   double tmp195_ = tmp189_ * tmp7_;
   double tmp196_ = tmp190_ * tmp40_;
   double tmp197_ = tmp195_ + tmp196_;
   double tmp198_ = (tmp197_) * (tmp29_);
   double tmp199_ = tmp194_ + tmp198_;
   double tmp200_ = tmp190_ * tmp2_;
   double tmp201_ = tmp200_ * (tmp37_);
   double tmp202_ = tmp199_ + tmp201_;
   double tmp203_ = tmp41_ * (tmp202_);
   double tmp204_ = -(tmp203_);
   double tmp205_ = tmp204_ / tmp86_;
   double tmp206_ = (tmp205_) * (tmp61_);
   double tmp207_ = tmp206_ * mLocPolyn2_State_0_0;
   double tmp208_ = (tmp207_) / tmp88_;
   double tmp209_ = (tmp208_) * (tmp66_);
   double tmp210_ = tmp190_ * tmp38_;
   double tmp211_ = tmp21_ * tmp39_;
   double tmp212_ = tmp210_ + tmp211_;
   double tmp213_ = (tmp212_) * (tmp17_);
   double tmp214_ = tmp190_ * tmp7_;
   double tmp215_ = tmp21_ * tmp40_;
   double tmp216_ = tmp214_ + tmp215_;
   double tmp217_ = (tmp216_) * (tmp29_);
   double tmp218_ = tmp213_ + tmp217_;
   double tmp219_ = tmp21_ * tmp2_;
   double tmp220_ = tmp219_ * (tmp37_);
   double tmp221_ = tmp218_ + tmp220_;
   double tmp222_ = (tmp221_) * (tmp54_);
   double tmp223_ = (tmp205_) * (tmp78_);
   double tmp224_ = tmp222_ + tmp223_;
   double tmp225_ = (tmp224_) * mLocPolyn2_State_0_0;
   double tmp226_ = (tmp225_) / tmp88_;
   double tmp227_ = (tmp226_) * (tmp83_);
   double tmp228_ = tmp108_ * (tmp44_);
   double tmp229_ = tmp41_ * tmp228_;
   double tmp230_ = -(tmp229_);
   double tmp231_ = tmp230_ / tmp86_;
   double tmp232_ = tmp108_ * tmp55_;
   double tmp233_ = tmp232_ * (tmp54_);
   double tmp234_ = (tmp231_) * (tmp61_);
   double tmp235_ = tmp233_ + tmp234_;
   double tmp236_ = (tmp235_) * mLocPolyn2_State_0_0;
   double tmp237_ = (tmp236_) / tmp88_;
   double tmp238_ = (tmp237_) * (tmp66_);
   double tmp239_ = tmp108_ * (tmp69_);
   double tmp240_ = tmp239_ * (tmp54_);
   double tmp241_ = (tmp231_) * (tmp78_);
   double tmp242_ = tmp240_ + tmp241_;
   double tmp243_ = (tmp242_) * mLocPolyn2_State_0_0;
   double tmp244_ = (tmp243_) / tmp88_;
   double tmp245_ = (tmp244_) * (tmp83_);
   double tmp246_ = tmp108_ * (tmp48_);
   double tmp247_ = tmp41_ * tmp246_;
   double tmp248_ = -(tmp247_);
   double tmp249_ = tmp248_ / tmp86_;
   double tmp250_ = tmp108_ * tmp57_;
   double tmp251_ = tmp250_ * (tmp54_);
   double tmp252_ = (tmp249_) * (tmp61_);
   double tmp253_ = tmp251_ + tmp252_;
   double tmp254_ = (tmp253_) * mLocPolyn2_State_0_0;
   double tmp255_ = (tmp254_) / tmp88_;
   double tmp256_ = (tmp255_) * (tmp66_);
   double tmp257_ = tmp108_ * (tmp73_);
   double tmp258_ = tmp257_ * (tmp54_);
   double tmp259_ = (tmp249_) * (tmp78_);
   double tmp260_ = tmp258_ + tmp259_;
   double tmp261_ = (tmp260_) * mLocPolyn2_State_0_0;
   double tmp262_ = (tmp261_) / tmp88_;
   double tmp263_ = (tmp262_) * (tmp83_);
   double tmp264_ = tmp108_ * tmp51_;
   double tmp265_ = tmp41_ * tmp264_;
   double tmp266_ = -(tmp265_);
   double tmp267_ = tmp266_ / tmp86_;
   double tmp268_ = tmp147_ * (tmp54_);
   double tmp269_ = (tmp267_) * (tmp61_);
   double tmp270_ = tmp268_ + tmp269_;
   double tmp271_ = (tmp270_) * mLocPolyn2_State_0_0;
   double tmp272_ = (tmp271_) / tmp88_;
   double tmp273_ = (tmp272_) * (tmp66_);
   double tmp274_ = tmp108_ * tmp76_;
   double tmp275_ = tmp274_ * (tmp54_);
   double tmp276_ = (tmp267_) * (tmp78_);
   double tmp277_ = tmp275_ + tmp276_;
   double tmp278_ = (tmp277_) * mLocPolyn2_State_0_0;
   double tmp279_ = (tmp278_) / tmp88_;
   double tmp280_ = (tmp279_) * (tmp83_);
   double tmp281_ = mLocProjI_x * (tmp44_);
   double tmp282_ = mLocProjI_y * (tmp48_);
   double tmp283_ = tmp281_ + tmp282_;
   double tmp284_ = mLocProjI_z * tmp51_;
   double tmp285_ = tmp283_ + tmp284_;
   double tmp286_ = tmp41_ * (tmp285_);
   double tmp287_ = -(tmp286_);
   double tmp288_ = tmp287_ / tmp86_;
   double tmp289_ = mLocProjI_x * tmp55_;
   double tmp290_ = mLocProjI_y * tmp57_;
   double tmp291_ = tmp289_ + tmp290_;
   double tmp292_ = mLocProjI_z * tmp8_;
   double tmp293_ = tmp291_ + tmp292_;
   double tmp294_ = (tmp293_) * (tmp54_);
   double tmp295_ = (tmp288_) * (tmp61_);
   double tmp296_ = tmp294_ + tmp295_;
   double tmp297_ = (tmp296_) * mLocPolyn2_State_0_0;
   double tmp298_ = (tmp297_) / tmp88_;
   double tmp299_ = (tmp298_) * (tmp66_);
   double tmp300_ = mLocProjI_x * (tmp69_);
   double tmp301_ = mLocProjI_y * (tmp73_);
   double tmp302_ = tmp300_ + tmp301_;
   double tmp303_ = mLocProjI_z * tmp76_;
   double tmp304_ = tmp302_ + tmp303_;
   double tmp305_ = (tmp304_) * (tmp54_);
   double tmp306_ = (tmp288_) * (tmp78_);
   double tmp307_ = tmp305_ + tmp306_;
   double tmp308_ = (tmp307_) * mLocPolyn2_State_0_0;
   double tmp309_ = (tmp308_) / tmp88_;
   double tmp310_ = (tmp309_) * (tmp83_);
   double tmp311_ = mLocProjJ_x * (tmp44_);
   double tmp312_ = mLocProjJ_y * (tmp48_);
   double tmp313_ = tmp311_ + tmp312_;
   double tmp314_ = mLocProjJ_z * tmp51_;
   double tmp315_ = tmp313_ + tmp314_;
   double tmp316_ = tmp41_ * (tmp315_);
   double tmp317_ = -(tmp316_);
   double tmp318_ = tmp317_ / tmp86_;
   double tmp319_ = mLocProjJ_x * tmp55_;
   double tmp320_ = mLocProjJ_y * tmp57_;
   double tmp321_ = tmp319_ + tmp320_;
   double tmp322_ = mLocProjJ_z * tmp8_;
   double tmp323_ = tmp321_ + tmp322_;
   double tmp324_ = (tmp323_) * (tmp54_);
   double tmp325_ = (tmp318_) * (tmp61_);
   double tmp326_ = tmp324_ + tmp325_;
   double tmp327_ = (tmp326_) * mLocPolyn2_State_0_0;
   double tmp328_ = (tmp327_) / tmp88_;
   double tmp329_ = (tmp328_) * (tmp66_);
   double tmp330_ = mLocProjJ_x * (tmp69_);
   double tmp331_ = mLocProjJ_y * (tmp73_);
   double tmp332_ = tmp330_ + tmp331_;
   double tmp333_ = mLocProjJ_z * tmp76_;
   double tmp334_ = tmp332_ + tmp333_;
   double tmp335_ = (tmp334_) * (tmp54_);
   double tmp336_ = (tmp318_) * (tmp78_);
   double tmp337_ = tmp335_ + tmp336_;
   double tmp338_ = (tmp337_) * mLocPolyn2_State_0_0;
   double tmp339_ = (tmp338_) / tmp88_;
   double tmp340_ = (tmp339_) * (tmp83_);
   double tmp341_ = mLocProjK_x * (tmp44_);
   double tmp342_ = mLocProjK_y * (tmp48_);
   double tmp343_ = tmp341_ + tmp342_;
   double tmp344_ = mLocProjK_z * tmp51_;
   double tmp345_ = tmp343_ + tmp344_;
   double tmp346_ = tmp41_ * (tmp345_);
   double tmp347_ = -(tmp346_);
   double tmp348_ = tmp347_ / tmp86_;
   double tmp349_ = mLocProjK_x * tmp55_;
   double tmp350_ = mLocProjK_y * tmp57_;
   double tmp351_ = tmp349_ + tmp350_;
   double tmp352_ = mLocProjK_z * tmp8_;
   double tmp353_ = tmp351_ + tmp352_;
   double tmp354_ = (tmp353_) * (tmp54_);
   double tmp355_ = (tmp348_) * (tmp61_);
   double tmp356_ = tmp354_ + tmp355_;
   double tmp357_ = (tmp356_) * mLocPolyn2_State_0_0;
   double tmp358_ = (tmp357_) / tmp88_;
   double tmp359_ = (tmp358_) * (tmp66_);
   double tmp360_ = mLocProjK_x * (tmp69_);
   double tmp361_ = mLocProjK_y * (tmp73_);
   double tmp362_ = tmp360_ + tmp361_;
   double tmp363_ = mLocProjK_z * tmp76_;
   double tmp364_ = tmp362_ + tmp363_;
   double tmp365_ = (tmp364_) * (tmp54_);
   double tmp366_ = (tmp348_) * (tmp78_);
   double tmp367_ = tmp365_ + tmp366_;
   double tmp368_ = (tmp367_) * mLocPolyn2_State_0_0;
   double tmp369_ = (tmp368_) / tmp88_;
   double tmp370_ = (tmp369_) * (tmp83_);
   double tmp371_ = 1 - tmp84_;
   double tmp372_ = (tmp92_) * (tmp83_);
   double tmp373_ = (tmp98_) * (tmp66_);
   double tmp374_ = tmp372_ + tmp373_;
   double tmp375_ = tmp100_ + tmp100_;
   double tmp376_ = tmp99_ * 2;
   double tmp377_ = tmp93_ + tmp93_;
   double tmp378_ = mCompCoord[8];
   double tmp379_ = (tmp102_) * tmp89_;
   double tmp380_ = tmp103_ + tmp103_;
   double tmp381_ = tmp104_ + tmp104_;
   double tmp382_ = (tmp66_) * mLocPolyn2_State_0_0;
   double tmp383_ = tmp382_ * mLocScNorm;
   double tmp384_ = tmp106_ * mLocPolyn2_State_0_0;
   double tmp385_ = tmp384_ * mLocScNorm;
   double tmp386_ = (tmp131_) * (tmp83_);
   double tmp387_ = (tmp145_) * (tmp66_);
   double tmp388_ = tmp386_ + tmp387_;
   double tmp389_ = tmp146_ + tmp146_;
   double tmp390_ = tmp132_ + tmp132_;
   double tmp391_ = (tmp173_) * (tmp83_);
   double tmp392_ = (tmp187_) * (tmp66_);
   double tmp393_ = tmp391_ + tmp392_;
   double tmp394_ = tmp188_ + tmp188_;
   double tmp395_ = tmp174_ + tmp174_;
   double tmp396_ = (tmp208_) * (tmp83_);
   double tmp397_ = (tmp226_) * (tmp66_);
   double tmp398_ = tmp396_ + tmp397_;
   double tmp399_ = tmp227_ + tmp227_;
   double tmp400_ = tmp209_ + tmp209_;
   double tmp401_ = (tmp237_) * (tmp83_);
   double tmp402_ = (tmp244_) * (tmp66_);
   double tmp403_ = tmp401_ + tmp402_;
   double tmp404_ = tmp245_ + tmp245_;
   double tmp405_ = tmp238_ + tmp238_;
   double tmp406_ = (tmp255_) * (tmp83_);
   double tmp407_ = (tmp262_) * (tmp66_);
   double tmp408_ = tmp406_ + tmp407_;
   double tmp409_ = tmp263_ + tmp263_;
   double tmp410_ = tmp256_ + tmp256_;
   double tmp411_ = (tmp272_) * (tmp83_);
   double tmp412_ = (tmp279_) * (tmp66_);
   double tmp413_ = tmp411_ + tmp412_;
   double tmp414_ = tmp280_ + tmp280_;
   double tmp415_ = tmp273_ + tmp273_;
   double tmp416_ = (tmp298_) * (tmp83_);
   double tmp417_ = (tmp309_) * (tmp66_);
   double tmp418_ = tmp416_ + tmp417_;
   double tmp419_ = tmp310_ + tmp310_;
   double tmp420_ = tmp299_ + tmp299_;
   double tmp421_ = (tmp328_) * (tmp83_);
   double tmp422_ = (tmp339_) * (tmp66_);
   double tmp423_ = tmp421_ + tmp422_;
   double tmp424_ = tmp340_ + tmp340_;
   double tmp425_ = tmp329_ + tmp329_;
   double tmp426_ = (tmp358_) * (tmp83_);
   double tmp427_ = (tmp369_) * (tmp66_);
   double tmp428_ = tmp426_ + tmp427_;
   double tmp429_ = tmp370_ + tmp370_;
   double tmp430_ = tmp359_ + tmp359_;

  mVal[0] = ((mLocPolyn2_State_1_0 + (((tmp85_) * (tmp66_) + tmp89_ * (tmp83_)) - tmp95_ * tmp105_ + tmp99_ * tmp106_ + tmp101_ * tmp107_) * mLocPolyn2_State_0_0) - mLocXIm) * mLocScNorm;

  mCompDer[0][0] = (((tmp92_) * (tmp85_) + (tmp98_) * tmp89_) - (tmp377_) * tmp95_ + (tmp374_) * tmp99_ + (tmp375_) * tmp101_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[0][1] = ((tmp102_) * (tmp85_) - (tmp380_) * tmp95_ + tmp104_ * tmp99_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[0][2] = (tmp379_ + tmp103_ * tmp99_ + (tmp381_) * tmp101_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[0][3] = tmp383_;
  mCompDer[0][4] = (tmp83_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[0][5] = -(2 * tmp105_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[0][6] = tmp385_;
  mCompDer[0][7] = tmp107_ * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[0][8] = 0;
  mCompDer[0][9] = (((tmp131_) * (tmp85_) + (tmp145_) * tmp89_) - (tmp390_) * tmp95_ + (tmp388_) * tmp99_ + (tmp389_) * tmp101_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[0][10] = (((tmp173_) * (tmp85_) + (tmp187_) * tmp89_) - (tmp395_) * tmp95_ + (tmp393_) * tmp99_ + (tmp394_) * tmp101_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[0][11] = (((tmp208_) * (tmp85_) + (tmp226_) * tmp89_) - (tmp400_) * tmp95_ + (tmp398_) * tmp99_ + (tmp399_) * tmp101_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[0][12] = (((tmp237_) * (tmp85_) + (tmp244_) * tmp89_) - (tmp405_) * tmp95_ + (tmp403_) * tmp99_ + (tmp404_) * tmp101_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[0][13] = (((tmp255_) * (tmp85_) + (tmp262_) * tmp89_) - (tmp410_) * tmp95_ + (tmp408_) * tmp99_ + (tmp409_) * tmp101_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[0][14] = (((tmp272_) * (tmp85_) + (tmp279_) * tmp89_) - (tmp415_) * tmp95_ + (tmp413_) * tmp99_ + (tmp414_) * tmp101_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[0][15] = (((tmp298_) * (tmp85_) + (tmp309_) * tmp89_) - (tmp420_) * tmp95_ + (tmp418_) * tmp99_ + (tmp419_) * tmp101_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[0][16] = (((tmp328_) * (tmp85_) + (tmp339_) * tmp89_) - (tmp425_) * tmp95_ + (tmp423_) * tmp99_ + (tmp424_) * tmp101_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[0][17] = (((tmp358_) * (tmp85_) + (tmp369_) * tmp89_) - (tmp430_) * tmp95_ + (tmp428_) * tmp99_ + (tmp429_) * tmp101_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mVal[1] = ((mLocPolyn2_State_2_0 + (((tmp371_) * (tmp83_) + tmp89_ * (tmp66_) + tmp94_ * tmp106_) - tmp376_ * tmp107_ + tmp378_ * tmp105_) * mLocPolyn2_State_0_0) - mLocYIm) * mLocScNorm;

  mCompDer[1][0] = (((tmp98_) * (tmp371_) + (tmp92_) * tmp89_ + (tmp374_) * tmp94_) - (tmp375_) * tmp376_ + (tmp377_) * tmp378_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[1][1] = (tmp379_ + tmp104_ * tmp94_ + (tmp380_) * tmp378_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[1][2] = (((tmp102_) * (tmp371_) + tmp103_ * tmp94_) - (tmp381_) * tmp376_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[1][3] = tmp108_ * (tmp83_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[1][4] = tmp383_;
  mCompDer[1][5] = tmp385_;
  mCompDer[1][6] = -(2 * tmp107_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[1][7] = 0;
  mCompDer[1][8] = tmp105_ * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[1][9] = (((tmp145_) * (tmp371_) + (tmp131_) * tmp89_ + (tmp388_) * tmp94_) - (tmp389_) * tmp376_ + (tmp390_) * tmp378_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[1][10] = (((tmp187_) * (tmp371_) + (tmp173_) * tmp89_ + (tmp393_) * tmp94_) - (tmp394_) * tmp376_ + (tmp395_) * tmp378_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[1][11] = (((tmp226_) * (tmp371_) + (tmp208_) * tmp89_ + (tmp398_) * tmp94_) - (tmp399_) * tmp376_ + (tmp400_) * tmp378_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[1][12] = (((tmp244_) * (tmp371_) + (tmp237_) * tmp89_ + (tmp403_) * tmp94_) - (tmp404_) * tmp376_ + (tmp405_) * tmp378_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[1][13] = (((tmp262_) * (tmp371_) + (tmp255_) * tmp89_ + (tmp408_) * tmp94_) - (tmp409_) * tmp376_ + (tmp410_) * tmp378_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[1][14] = (((tmp279_) * (tmp371_) + (tmp272_) * tmp89_ + (tmp413_) * tmp94_) - (tmp414_) * tmp376_ + (tmp415_) * tmp378_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[1][15] = (((tmp309_) * (tmp371_) + (tmp298_) * tmp89_ + (tmp418_) * tmp94_) - (tmp419_) * tmp376_ + (tmp420_) * tmp378_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[1][16] = (((tmp339_) * (tmp371_) + (tmp328_) * tmp89_ + (tmp423_) * tmp94_) - (tmp424_) * tmp376_ + (tmp425_) * tmp378_) * mLocPolyn2_State_0_0 * mLocScNorm;
  mCompDer[1][17] = (((tmp369_) * (tmp371_) + (tmp358_) * tmp89_ + (tmp428_) * tmp94_) - (tmp429_) * tmp376_ + (tmp430_) * tmp378_) * mLocPolyn2_State_0_0 * mLocScNorm;
}


void cEqAppui_PProjInc_M2CPolyn2::ComputeValDerivHessian()
{
  ELISE_ASSERT(false,"Foncteur cEqAppui_PProjInc_M2CPolyn2 Has no Der Sec");
}

void cEqAppui_PProjInc_M2CPolyn2::SetPolyn2_State_0_0(double aVal){ mLocPolyn2_State_0_0 = aVal;}
void cEqAppui_PProjInc_M2CPolyn2::SetPolyn2_State_1_0(double aVal){ mLocPolyn2_State_1_0 = aVal;}
void cEqAppui_PProjInc_M2CPolyn2::SetPolyn2_State_2_0(double aVal){ mLocPolyn2_State_2_0 = aVal;}
void cEqAppui_PProjInc_M2CPolyn2::SetProjI_x(double aVal){ mLocProjI_x = aVal;}
void cEqAppui_PProjInc_M2CPolyn2::SetProjI_y(double aVal){ mLocProjI_y = aVal;}
void cEqAppui_PProjInc_M2CPolyn2::SetProjI_z(double aVal){ mLocProjI_z = aVal;}
void cEqAppui_PProjInc_M2CPolyn2::SetProjJ_x(double aVal){ mLocProjJ_x = aVal;}
void cEqAppui_PProjInc_M2CPolyn2::SetProjJ_y(double aVal){ mLocProjJ_y = aVal;}
void cEqAppui_PProjInc_M2CPolyn2::SetProjJ_z(double aVal){ mLocProjJ_z = aVal;}
void cEqAppui_PProjInc_M2CPolyn2::SetProjK_x(double aVal){ mLocProjK_x = aVal;}
void cEqAppui_PProjInc_M2CPolyn2::SetProjK_y(double aVal){ mLocProjK_y = aVal;}
void cEqAppui_PProjInc_M2CPolyn2::SetProjK_z(double aVal){ mLocProjK_z = aVal;}
void cEqAppui_PProjInc_M2CPolyn2::SetProjP0_x(double aVal){ mLocProjP0_x = aVal;}
void cEqAppui_PProjInc_M2CPolyn2::SetProjP0_y(double aVal){ mLocProjP0_y = aVal;}
void cEqAppui_PProjInc_M2CPolyn2::SetProjP0_z(double aVal){ mLocProjP0_z = aVal;}
void cEqAppui_PProjInc_M2CPolyn2::SetScNorm(double aVal){ mLocScNorm = aVal;}
void cEqAppui_PProjInc_M2CPolyn2::SetXIm(double aVal){ mLocXIm = aVal;}
void cEqAppui_PProjInc_M2CPolyn2::SetYIm(double aVal){ mLocYIm = aVal;}



double * cEqAppui_PProjInc_M2CPolyn2::AdrVarLocFromString(const std::string & aName)
{
   if (aName == "Polyn2_State_0_0") return & mLocPolyn2_State_0_0;
   if (aName == "Polyn2_State_1_0") return & mLocPolyn2_State_1_0;
   if (aName == "Polyn2_State_2_0") return & mLocPolyn2_State_2_0;
   if (aName == "ProjI_x") return & mLocProjI_x;
   if (aName == "ProjI_y") return & mLocProjI_y;
   if (aName == "ProjI_z") return & mLocProjI_z;
   if (aName == "ProjJ_x") return & mLocProjJ_x;
   if (aName == "ProjJ_y") return & mLocProjJ_y;
   if (aName == "ProjJ_z") return & mLocProjJ_z;
   if (aName == "ProjK_x") return & mLocProjK_x;
   if (aName == "ProjK_y") return & mLocProjK_y;
   if (aName == "ProjK_z") return & mLocProjK_z;
   if (aName == "ProjP0_x") return & mLocProjP0_x;
   if (aName == "ProjP0_y") return & mLocProjP0_y;
   if (aName == "ProjP0_z") return & mLocProjP0_z;
   if (aName == "ScNorm") return & mLocScNorm;
   if (aName == "XIm") return & mLocXIm;
   if (aName == "YIm") return & mLocYIm;
   return 0;
}


cElCompiledFonc::cAutoAddEntry cEqAppui_PProjInc_M2CPolyn2::mTheAuto("cEqAppui_PProjInc_M2CPolyn2",cEqAppui_PProjInc_M2CPolyn2::Alloc);


cElCompiledFonc *  cEqAppui_PProjInc_M2CPolyn2::Alloc()
{  return new cEqAppui_PProjInc_M2CPolyn2();
}

