
/**
 * Basic constants of VCF files
 */
public class VCFconstants {
	/**
	 * SV types: represented in our genome graph.
	 *
	 * Remark: for simplicity, an INVDUP is assumed to be a sequence $xVV'y$, 
	 * where $V'$ is the reverse-complement of $V$ and the reference is $xVy$. 
	 * This is not equivalent to Sniffles' definition (Supplementary Figure 
	 * 2.4), in which it seems to be $xVy'V'z$ where the reference is $xVyz$. So 
	 * we are interpreting Sniffles' INVDUP as if it were Sniffles' DUP+INV 
	 * (see again Supplementary Figure 2.4).
	 *
	 * Remark: a DEL/INV is assumed to be a sequence $aC'e$ where $C'$ is the
	 * reverse-complement of $C$ and the reference is $aBCDe$. The length of B,D
	 * is unknown, so we add to the graph every possible edge between a position
	 * in the first half and a position in the second half of BCD.
	 */
	public static final byte TYPE_INSERTION = 1;
	public static final byte TYPE_DELETION = 2;
	public static final byte TYPE_DEL_INV = 3;
	public static final byte TYPE_INVERSION = 4;
	public static final byte TYPE_INV_DUP = 5;
	public static final byte TYPE_DUPLICATION = 6;
	public static final byte TYPE_CNV = 7;
	public static final byte TYPE_BREAKEND = 8;
	public static final byte TYPE_TRANSLOCATION = 9;
	
	/**
	 * SV types: labels used by callers.
	 */
	public static final String DEL_STR = "DEL";
	public static final String DEL_ME_STR = "DEL:ME";
	public static final String DEL_INV_STR = "DEL/INV";
	public static final String INS_STR = "INS";
	public static final String INS_ME_STR = "INS:ME";
	public static final String INS_NOVEL_STR = "INS:NOVEL";
	public static final String DUP_STR = "DUP";
	public static final String DUP_TANDEM_STR = "DUP:TANDEM";
	public static final String DUP_INT_STR = "DUP:INT";
	public static final String INV_STR = "INV";
	public static final String INV_DUP_STR = "INVDUP";
	public static final String CNV_STR = "CNV";
	public static final String BND_STR = "BND";
	public static final String TRA_STR = "TRA";
	
	/**
	 * Basic constants
	 */
	public static final char COMMENT = '#';
	public static final String SEPARATOR = ";";
	public static final String PASS_STR = "PASS";
	public static final String PRECISE_STR = "PRECISE";
	public static final String IMPRECISE_STR = "IMPRECISE";
	public static final String END_STR = "END";
	public static final String SVTYPE_STR = "SVTYPE";
	public static final String SVLEN_STR = "SVLEN";
	public static final String CHR2_STR = "CHR2";
	public static final String CT_STR = "CT";
	public static final String CT_325_STR = "'3to5'";
	public static final String CT_523_STR = "'5to3'";
	public static final String CT_525_STR = "'5to5'";
	public static final String CT_323_STR = "'3to3'";
	
	/**
	 * Confidence intervals of positions.
	 *
	 * Remark: some callers report a standard deviation instead of a confidence
	 * interval. Sniffles reports additional interval information in its BEDPE
	 * output (an alternative to VCF), which our programs disregard for 
	 * simplicity. Some callers use CILEN to express a "confidence interval 
	 * around inserted/deleted material between breakends": we interpret CILEN
	 * exactly like CIEND, and we ignore it for insertions (since representing
	 * variable-length insertion strings complicates our code).
	 */
	public static final String CI_SEPARATOR = ",";
	public static final String CIPOS_STR = "CIPOS";
	public static final String CIEND_STR = "CIEND";
	public static final String STD_START1_STR = "STD_quant_start";
	public static final String STD_START2_STR = "STD_POS1";
	public static final String STD_END1_STR = "STD_quant_stop";
	public static final String STD_END2_STR = "STD_POS2";
	public static final String CILEN_STR = "CILEN";
	public static final int CIPOS_STR_LENGTH = CIPOS_STR.length();
	public static final int CIEND_STR_LENGTH = CIEND_STR.length();
	public static final int STD_START1_STR_LENGTH = STD_START1_STR.length();
	public static final int STD_START2_STR_LENGTH = STD_START2_STR.length();
	public static final int STD_END1_STR_LENGTH = STD_END1_STR.length();
	public static final int STD_END2_STR_LENGTH = STD_END2_STR.length();
	public static final int CILEN_STR_LENGTH = CILEN_STR.length();
	
	/**
	 * Chromosomes
	 */
	public static final int N_CHROMOSOMES = 25;  // 22+X+Y+M
	public static final String CHR_STR = "chr";
	public static final String X_STR_PRIME = "X";
	public static final String Y_STR_PRIME = "Y";
	public static final String M_STR_PRIME = "M";
	public static final String MT_STR_PRIME = "MT";
	public static final int CHR_STR_LENGTH = CHR_STR.length();
	public static final String X_STR = CHR_STR+X_STR_PRIME;
	public static final String Y_STR = CHR_STR+Y_STR_PRIME;
	public static final String M_STR = CHR_STR+M_STR_PRIME;
	public static final String MT_STR = CHR_STR+MT_STR_PRIME;
	
	
	/**
	 * @return NULL if $field$ does not occur in $str$.
	 */
	public static final String getField(String str, String field) {
		final int FIELD_LENGTH = field.length()+1;
		int p = str.indexOf(field+"=");
		if (p<0) return null;
		if (field.equalsIgnoreCase(END_STR)) {
			while (p>=2 && str.substring(p-2,p-2+CIEND_STR.length()).equalsIgnoreCase(CIEND_STR)) p=str.indexOf(field+"=",p+1);
			if (p<0) return null;
		}
		final int q = str.indexOf(SEPARATOR,p+FIELD_LENGTH);
		return str.substring(p+FIELD_LENGTH,q<0?str.length():q);
	}
	
	
	/**
	 * @return -1 iff the type cannot be determined.
	 */
	public static final int getType_infoField(String type) {
		if (type==null || type.length()==0) return -1;
		if ( type.equalsIgnoreCase(DEL_STR) || 
			 type.equalsIgnoreCase(DEL_ME_STR)
		   ) return TYPE_DELETION;
		else if (type.equalsIgnoreCase(DEL_INV_STR)) return TYPE_DEL_INV;
		else if ( type.equalsIgnoreCase(INS_STR) || 
			      type.equalsIgnoreCase(INS_ME_STR) || 
				  type.equalsIgnoreCase(INS_NOVEL_STR)
				) return TYPE_INSERTION;
		else if ( type.equalsIgnoreCase(DUP_STR) ||
			      type.equalsIgnoreCase(DUP_TANDEM_STR) ||
				  type.equalsIgnoreCase(DUP_INT_STR)
			    ) return TYPE_DUPLICATION;
		else if (type.equalsIgnoreCase(INV_STR)) return TYPE_INVERSION;
		else if (type.equalsIgnoreCase(INV_DUP_STR)) return TYPE_INV_DUP;
		else if (type.equalsIgnoreCase(CNV_STR)) return TYPE_CNV;
		else if (type.equalsIgnoreCase(BND_STR)) return TYPE_BREAKEND;
		else if (type.equalsIgnoreCase(TRA_STR)) return TYPE_TRANSLOCATION;
		else return -1;
	}
	
	
	public static final int getType_altField(String str) {
		if (str==null || str.length()==0) return -1;
		if ( str.equalsIgnoreCase("<"+DEL_STR+">") || 
			 str.equalsIgnoreCase("<"+DEL_ME_STR+">")
		   ) return TYPE_DELETION;
		else if (str.equalsIgnoreCase("<"+DEL_INV_STR+">")) return TYPE_DEL_INV;
		else if ( str.equalsIgnoreCase("<"+INS_STR+">") ||
			      str.equalsIgnoreCase("<"+INS_ME_STR+">") ||
				  str.equalsIgnoreCase("<"+INS_NOVEL_STR+">")
			    ) return TYPE_INSERTION;
		else if ( str.equalsIgnoreCase("<"+DUP_STR+">") ||
			      str.equalsIgnoreCase("<"+DUP_TANDEM_STR+">") ||
				  str.equalsIgnoreCase("<"+DUP_INT_STR+">")
			    ) return TYPE_DUPLICATION;
		else if (str.equalsIgnoreCase("<"+INV_STR+">")) return TYPE_INVERSION;
		else if (str.equalsIgnoreCase("<"+INV_DUP_STR+">")) return TYPE_INV_DUP;
		else if (str.equalsIgnoreCase("<"+CNV_STR+">")) return TYPE_CNV;
		else if (str.equalsIgnoreCase("<"+BND_STR+">")) return TYPE_BREAKEND;
		else if (str.equalsIgnoreCase("<"+TRA_STR+">")) return TYPE_TRANSLOCATION;
		else return -1;
	}
	
	
	/**
	 * @return TRUE iff $str$ is a supported SV type.
	 */
	public static final boolean isType(String str) {
		return str.equalsIgnoreCase(DEL_STR) || str.equalsIgnoreCase(DEL_ME_STR) || 
			   str.equalsIgnoreCase(DEL_INV_STR) ||
			   str.equalsIgnoreCase(INS_STR) || str.equalsIgnoreCase(INS_ME_STR) || str.equalsIgnoreCase(INS_NOVEL_STR) ||
			   str.equalsIgnoreCase(DUP_STR) || str.equalsIgnoreCase(DUP_TANDEM_STR) || str.equalsIgnoreCase(DUP_INT_STR) ||
			   str.equalsIgnoreCase(INV_STR) || 
			   str.equalsIgnoreCase(INV_DUP_STR) ||
			   str.equalsIgnoreCase(CNV_STR) ||
			   str.equalsIgnoreCase(BND_STR) ||
			   str.equalsIgnoreCase(TRA_STR);
	}
	
	
	/**
	 * @param out output array: 0=left side, 1=right side.
	 */
	public static final void ct2side(String str, int[] out) {
		if (str.equalsIgnoreCase(CT_325_STR)) { out[0]=1; out[1]=0; }
		else if (str.equalsIgnoreCase(CT_523_STR)) { out[0]=0; out[1]=1; }
		else if (str.equalsIgnoreCase(CT_525_STR)) { out[0]=0; out[1]=0; }
		else if (str.equalsIgnoreCase(CT_323_STR)) { out[0]=1; out[1]=1; }	
	} 
	
	
	/**
	 * @return TRUE iff $str$ represents the ID of a chromosome.
	 */
	public static final boolean isChromosome(String str) {
		if ( Character.isDigit(str.charAt(0)) ||
			 str.equalsIgnoreCase(X_STR) || str.equalsIgnoreCase(X_STR_PRIME) || 
			 str.equalsIgnoreCase(Y_STR) || str.equalsIgnoreCase(Y_STR_PRIME) ||
			 str.equalsIgnoreCase(M_STR) || str.equalsIgnoreCase(M_STR_PRIME) ||
			 str.equalsIgnoreCase(MT_STR) || str.equalsIgnoreCase(MT_STR_PRIME) ||
			 (str.substring(0,CHR_STR_LENGTH).equalsIgnoreCase(CHR_STR) && str.length()<=CHR_STR_LENGTH+2)
		   ) return true;
		return false;
	}
	
	
	/**
	 * Simplified variant of $VCF2gfa.string2contig()$ that only works for
	 * chromosomes.
	 *
	 * @return one-based.
	 */
	public static final int string2contig(String str) {
		if (str.length()>=CHR_STR_LENGTH && str.substring(0,CHR_STR_LENGTH).equalsIgnoreCase(CHR_STR)) {
			if (str.equalsIgnoreCase(X_STR)) return 23;
			else if (str.equalsIgnoreCase(Y_STR)) return 24;
			else if (str.equalsIgnoreCase(M_STR) || str.equalsIgnoreCase(MT_STR)) return 25;
			else if (str.length()<=CHR_STR_LENGTH+2) return Integer.parseInt(str.substring(CHR_STR_LENGTH));
			else return -1;
		}
		else {
			if (str.equalsIgnoreCase(X_STR_PRIME)) return 23;
			else if (str.equalsIgnoreCase(Y_STR_PRIME)) return 24;
			else if (str.equalsIgnoreCase(M_STR_PRIME) || str.equalsIgnoreCase(MT_STR_PRIME)) return 25;
			else if (str.length()<=2) return Integer.parseInt(str);
			else return -1;
		}
	}
	
}