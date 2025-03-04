import sys
import copy

indels = set()
asv_set = set()

with open(sys.argv[1],'r') as asv:
    for l in asv:
        l = l.strip().split(',')
        indels.add(l[1])
        asv_set.add(l[0])

indel_dict = {}
for i in indels:
    indel_dict[i] = str(0)

asv_dict = {}
for a in asv_set:
    asv_dict[a] = copy.deepcopy(indel_dict)

with open(sys.argv[1],'r') as asv:
    for l in asv:
        l = l.strip().split(',')
        indel = l[1]
        asv_name = l[0]
        asv_dict[asv_name][indel] = str(1)
header = "ASV" + ',' + ','.join(list(indel_dict.keys()))
print(header)
for k,v in asv_dict.items():
    outline = k + ',' + ','.join(list(v.values()))
    print(outline)


#ASV694,d_147_30nts
#ASV694,d_17_1nts
#ASV695,d_146_2nts
#ASV695,d_17_1nts


#ASV,d_17_1nts,d_145_8nts,i_147_1nts,d_147_1nts,i_20_5nts,i_44_4nts,d_16_2nts,i_147_2nts,d_146_2nts,d_147_5nts,d_17_14nts,i_43_2nts,i_149_1nts,i_146_1nts,d_19_1nts,d_144_4nts,d_12_9nts,d_143_5nts,d_148_5nts,d_42_2nts,d_46_4nts,d_17_11nts,d_19_10nts,d_140_26nts,d_18_2nts,d_148_4nts,d_173_4nts,d_44_2nts,d_150_1nts,i_146_8nts,d_17_2nts,d_139_11nts,d_148_25nts,d_9_12nts,d_13_8nts,d_144_40nts,d_146_3nts,d_40_4nts,d_143_6nts,d_172_2nts,d_145_6nts,d_141_6nts,i_19_3nts,d_174_17nts,d_173_3nts,d_16_3nts,d_12_7nts,d_94_3nts,d_13_14nts,d_136_39nts,d_18_9nts,d_7_19nts,d_8_26nts,d_8_50nts,d_18_25nts,d_103_130nts,d_137_62nts,d_177_52nts,d_149_1nts,d_164_52nts,d_147_10nts,d_147_30nts,d_138_50nts,d_148_10nts,d_149_7nts,d_164_12nts,d_137_58nts,d_45_1nts,d_148_11nts,d_162_14nts,d_128_52nts,d_190_21nts,d_144_49nts,d_148_24nts,d_150_26nts,d_149_27nts,d_147_4nts,d_146_29nts,d_146_37nts,d_45_3nts,d_151_26nts,d_141_36nts,d_43_1nts,d_11_10nts,d_8_78nts,i_95_1nts,d_8_28nts,d_8_102nts,d_8_156nts,d_8_10nts,d_18_129nts,d_8_130nts,d_8_11nts,d_142_31nts,d_9_15nts,d_134_30nts,d_20_25nts,d_34_15nts,d_10_17nts,d_10_54nts,d_144_34nts,d_146_28nts,d_11_6nts,d_10_210nts,d_133_40nts,d_148_3nts,d_145_2nts,d_149_10nts,d_166_10nts,i_16_7nts,i_20_2nts,d_16_4nts,d_12_6nts,d_174_11nts,d_139_63nts,d_13_4nts,d_148_38nts,d_147_16nts,d_150_25nts,i_146_3nts,d_67_8nts,d_18_46nts,d_66_28nts,d_13_10nts,d_141_12nts,d_14_4nts,d_14_2nts,d_14_16nts,d_51_130nts,d_14_164nts,d_17_5nts,d_18_5nts,d_18_7nts,d_19_213nts,d_18_24nts,i_17_2nts,d_18_128nts,d_148_2nts,d_148_29nts,d_18_8nts,d_129_22nts,d_153_55nts,i_146_6nts,d_135_50nts,d_19_2nts,d_20_7nts,i_17_1nts,d_145_5nts,d_147_28nts,d_128_25nts,d_102_78nts,d_95_8nts,d_18_28nts,i_43_1nts,d_50_26nts,d_148_14nts,i_19_4nts,d_18_3nts,d_72_130nts,i_18_2nts,d_20_26nts,d_20_130nts,d_80_52nts,d_141_18nts,d_161_16nts,d_137_14nts,d_46_130nts,d_36_26nts,d_30_176nts,d_24_156nts,d_19_29nts,i_18_30nts,i_19_9nts,d_18_4nts,d_18_133nts,d_18_30nts,i_17_20nts,i_17_113nts,d_17_18nts,d_17_27nts,d_17_30nts,i_18_5nts,d_15_5nts,d_15_37nts,d_12_195nts,d_122_90nts,d_11_38nts,i_14_5nts,d_9_40nts,d_9_30nts,d_6_13nts,d_39_10nts,d_42_5nts,d_93_19nts,d_143_7nts,i_173_1nts,i_19_2nts,d_124_26nts,d_144_30nts,d_136_47nts,d_129_20nts,d_145_27nts,d_44_4nts,d_118_37nts,d_43_4nts,d_140_10nts,i_146_10nts,d_147_40nts,d_8_37nts,d_137_20nts,d_105_52nts,d_144_3nts,d_148_44nts,d_134_23nts,d_141_10nts,d_138_15nts,d_108_57nts,d_18_18nts,d_68_17nts,d_91_11nts,d_19_5nts,d_19_39nts,i_19_10nts,d_140_7nts,d_144_11nts,d_19_6nts,d_17_16nts,d_17_79nts,d_17_56nts,d_17_134nts,d_60_14nts,d_16_7nts,d_18_6nts,i_19_18nts,d_42_3nts,d_15_11nts,d_39_4nts,d_16_6nts,i_121_4nts,d_45_15nts,d_125_78nts,d_153_27nts,d_46_122nts,d_16_132nts,d_132_26nts,d_10_80nts,d_162_23nts,d_97_1nts,d_13_6nts,d_13_140nts,d_129_18nts,d_15_19nts,d_142_20nts,d_143_36nts,d_18_12nts,d_18_14nts,d_48_135nts,d_146_44nts,d_150_15nts,i_69_4nts,d_17_189nts,d_16_28nts,d_16_13nts,d_11_12nts,d_10_19nts
#ASV556,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
