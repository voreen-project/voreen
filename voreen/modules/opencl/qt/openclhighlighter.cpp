/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

// This class is based on the OpenCL syntax highlighter written by Markus Kramer for
// Shader Maker - a cross-platform OpenCL editor

#include "modules/opencl/qt/openclhighlighter.h"

namespace voreen {

OpenCLHighlighter::OpenCLHighlighter(QTextDocument* doc)
    : SyntaxHighlighter(doc)
{
    setupKeywords();
    setupReservedKeywords();
    setupPreprocessor();
    setupNumberRules();
    setupRegisters();
    setupBuiltInFunctions();
    setupComments();
}

void OpenCLHighlighter::setupComments() {
    highlightRule_t rule;

    m_singleLineCommentFormat.setForeground(Qt::darkGreen);
    m_singleLineCommentFormat.setFontItalic(true);

    rule.pattern = QRegExp("//[^\n]*");
    rule.format = m_singleLineCommentFormat;
    m_rules.append(rule);

    m_multiLineCommentFormat.setForeground(Qt::darkGreen);
    m_multiLineCommentFormat.setFontItalic(true);

    m_commentStartExpression = QRegExp("/\\*");
    m_commentEndExpression = QRegExp("\\*/");
}

void OpenCLHighlighter::setupPreprocessor() {
    m_preprocessorFormat.setForeground(Qt::blue);

    // highlight each line beginning with #
    highlightRule_t rule;
    rule.pattern = QRegExp("#[^\n]*");
    rule.format = m_preprocessorFormat;
    m_rules.append(rule);
}

void OpenCLHighlighter::setupNumberRules() {
    highlightRule_t rule;

    m_integerFormat.setForeground(Qt::magenta);
    m_floatFormat.setForeground(Qt::darkBlue);

    rule.pattern = QRegExp("\\b[0-9]+\\b");
    rule.format = m_integerFormat;
    m_rules.append(rule);

    // must be appended after integer rule to "overwrite" it.
    rule.pattern = QRegExp("\\b[0-9]\\.[0-9]+\\b");
    rule.format = m_floatFormat;
    m_rules.append(rule);
}

void OpenCLHighlighter::setupRegisters() {
    highlightRule_t rule;
    QStringList keywordPatterns;

    m_registerFormat.setForeground(Qt::darkBlue);
    m_registerFormat.setFontWeight(QFont::Bold);

    //
    // casts and constness
    //
    keywordPatterns
        << "\\bvec_type_hint\\b" << "\\bwork_group_size_hint\\b" << "\\baligned packed endian\\b"
        << "\\bas_(uchar|char|ushort|short|uint|int|ulong|long|float|double)\\b"
        << "\\bas_(uchar|char|ushort|short|uint|int|ulong|long|float|double)(2|3|4|8|16)\\b"
        << "\\bconvert_(uchar|char|ushort|short|uint|int|ulong|long|float|double)(2|3|4|8|16)?\\b"
        << "\\bconvert_(uchar|char|ushort|short|uint|int|ulong|long|float|double)(2|3|4|8|16)?_sat\\b"
        << "\\bconvert_(uchar|char|ushort|short|uint|int|ulong|long|float|double)(2|3|4|8|16)?_sat_(rte|rtz|rtp|rtn)\\b"
        << "\\bglobal\\b" << "\\b__global\\b" << "\\blocal\\b" << "\\b__local\\b" << "\\bconstant\\b" << "\\b__constant\\b" << "\\bprivate\\b" << "\\b__private\\b"
        << "\\bkernel\\b" << "\\b__kernel\\b" << "\\b __attribute__\\b"
        << "\\bread_only\\b" << "\\b__read_only\\b" << "\\bwrite_only\\b" << "\\b__write_only\\b"
        << "\\bcomplex\\b" << "\\bimaginary\\b"
        ;

    foreach(QString pattern, keywordPatterns) {
        rule.pattern = QRegExp(pattern);
        rule.format = m_registerFormat;
        m_rules.append(rule);
    }
}

void OpenCLHighlighter::setupBuiltInFunctions() {
    highlightRule_t rule;
    QStringList keywordPatterns;

    m_builtInFunctionFormat.setForeground(Qt::darkBlue);
    m_builtInFunctionFormat.setFontWeight(QFont::Bold);

    keywordPatterns // trigonometry
        << "\\bget_work_dim\\b" << "\\bget_global_size\\b"
        << "\\bget_global_id\\b" << "\\bget_local_size\\b" << "\\bget_local_id\\b"
        << "\\bget_num_groups\\b" << "\\bget_group_id\\b"
        << "\\bget_global_offset\\b" << "\\bcos\\b" << "\\bcosh\\b"
        << "\\bcospi\\b" << "\\bacos\\b" << "\\bacosh\\b" << "\\bacospi\\b"
        << "\\bsin\\b" << "\\bsincos\\b" << "\\bsinh\\b" << "\\bsinpi\\b"
        << "\\basin\\b" << "\\basinh\\b" << "\\basinpi\\b" << "\\btan\\b"
        << "\\btanh\\b" << "\\btanpi\\b" << "\\batan\\b" << "\\batan2\\b"
        << "\\batanh\\b" << "\\batanpi\\b" << "\\batan2pi\\b" << "\\bcbrt\\b"
        << "\\bceil\\b" << "\\bcopysign\\b" << "\\berfc\\b" << "\\berf\\b"
        << "\\bexp\\b" << "\\bexp2\\b" << "\\bexp10\\b" << "\\bexpm1\\b"
        << "\\bfabs\\b" << "\\bfdim\\b" << "\\bfloor\\b" << "\\bfma\\b"
        << "\\bfmax\\b" << "\\bfmin\\b" << "\\bfract\\b" << "\\bfrexp\\b"
        << "\\bhypot\\b" << "\\bilogb\\b" << "\\bldexp\\b" << "\\bldexp\\b"
        << "\\blgamma\\b" << "\\blgamma_r\\b" << "\\blog\\b" << "\\blog2\\b"
        << "\\blog10\\b" << "\\blog1p\\b" << "\\blogb\\b" << "\\bmad\\b"
        << "\\bmodf\\b" << "\\bnan\\b" << "\\bnextafter\\b" << "\\bpow\\b"
        << "\\bpown\\b" << "\\bpowr\\b" << "\\bremainder\\b" << "\\bremquo\\b"
        << "\\brint\\b" << "\\brootn\\b" << "\\bround\\b" << "\\brsqrt\\b"
        << "\\bsqrt\\b" << "\\btgamma\\b" << "\\btrunc\\b" << "\\bhalf_cos\\b"
        << "\\bhalf_divide\\b" << "\\bhalf_exp\\b" << "\\bhalf_exp2\\b"
        << "\\bhalf_exp10\\b" << "\\bhalf_log\\b" << "\\bhalf_log2\\b"
        << "\\bhalf_log10\\b" << "\\bhalf_powr\\b" << "\\bhalf_recip\\b"
        << "\\bhalf_rsqrt\\b" << "\\bhalf_sin\\b" << "\\bhalf_sqrt\\b"
        << "\\bhalf_tan\\b" << "\\bnative_cos\\b" << "\\bnative_divide\\b"
        << "\\bnative_exp\\b" << "\\bnative_exp2\\b" << "\\bnative_exp10\\b"
        << "\\bnative_log\\b" << "\\bnative_log2\\b" << "\\bnative_log10\\b"
        << "\\bnative_powr\\b" << "\\bnative_recip\\b" << "\\bnative_rsqrt\\b"
        << "\\bnative_sin\\b" << "\\bnative_sqrt\\b" << "\\bnative_tan\\b"
        << "\\babs\\b" << "\\babs_diff\\b" << "\\badd_sat\\b" << "\\bhadd\\b"
        << "\\brhadd\\b" << "\\bclz\\b" << "\\bmad_hi\\b" << "\\bmad_sat\\b"
        << "\\bmax\\b" << "\\bmin\\b" << "\\bmul_hi\\b" << "\\brotate\\b"
        << "\\bsub_sat\\b" << "\\bupsample\\b" << "\\bmad24\\b" << "\\bmul24\\b"
        << "\\bclamp\\b" << "\\bdegrees\\b" << "\\bmax\\b" << "\\bmin\\b"
        << "\\bmix\\b" << "\\bradians\\b" << "\\bstep\\b" << "\\bsmoothstep\\b"
        << "\\bsign\\b" << "\\bcross\\b" << "\\bdot\\b" << "\\bdistance\\b"
        << "\\blength\\b" << "\\bnormalize\\b" << "\\bfast_distance\\b"
        << "\\bfast_length\\b" << "\\bfast_normalize\\b" << "\\bvec_step\\b"
        << "\\bshuffle\\b" << "\\bshuffle2\\b" << "\\bisequal\\b"
        << "\\bisnotequal\\b" << "\\bisgreater\\b" << "\\bisgreaterequal\\b"
        << "\\bisless\\b" << "\\bislessequal\\b" << "\\bislessgreater\\b"
        << "\\bisfinite\\b" << "\\bisinf\\b" << "\\bisnan\\b" << "\\bisnormal\\b"
        << "\\bisordered\\b" << "\\bisunordered\\b" << "\\bsignbit\\b"
        << "\\bany\\b" << "\\ball\\b" << "\\bbitselect\\b" << "\\bselect\\b"
        << "\\bvload_half\\b" << "\\bvstore_half\\b" << "\\bget_image_width\\b"
        << "\\bget_image_height\\b" << "\\bget_image_depth\\b"
        << "\\bget_image_channel_data_type\\b" << "\\bget_image_channel_order\\b"
        << "\\bget_image_dim\\b" << "\\bbarrier\\b" << "\\bmem_fence\\b"
        << "\\bread_mem_fence\\b" << "\\bwrite_mem_fence\\b"
        << "\\ba<<c_work_group_copy\\b" << "\\ba<<c_work_group__strided_copy\\b"
        << "\\bwait_group_events\\b" << "\\bprefetch\\b"
        << "\\bvload(2|3|4|8|16)\\b"
        << "\\bvload_half(2|3|4|8|16)?\\b"
        << "\\bvloada_half(2|3|4|8|16)\\b"
        << "\\bvloada_half(2|3|4|8|16)_(rte|rtz|rtp)\\b"
        << "\\bvstore(2|3|4|8|16)\\b"
        << "\\bvstore(rte|rtz|rtp|rtn)\\b"
        << "\\bvstore_half(2|3|4|8|16)?\\b"
        << "\\bvstore_half_(rte|rtz|rtp|rtn)\\b"
        << "\\bvstore_half(2|3|4|8|16_(rte|rtz|rtp|rtn)\\b"
        << "\\bvstorea_half(2|3|4|8|16)\\b"
        << "\\bvstorea_half(2|3|4|8|16)_(rte|rtz|rtp|rtn)\\b"
        << "\\bread_image(f|i|ui|h)\\b"
        << "\\bwrite_image(f|i|ui|h)\\b"
        << "\\batom_(add|sub|xchg|inc|dec|cmpxchg|min|max|and|or|xor)\\b"
        ;

    foreach(QString pattern, keywordPatterns) {
        rule.pattern = QRegExp(pattern);
        rule.format = m_builtInFunctionFormat;
        m_rules.append(rule);
    }
}

void OpenCLHighlighter::setupKeywords() {
    highlightRule_t rule;
    QStringList keywordPatterns;

    m_keywordFormat.setForeground(Qt::blue);

    //
    // used keywords
    //
    keywordPatterns << "\\bconst\\b" << "\\bbreak\\b" << "\\bcontinue\\b"
                    << "\\bdo\\b" << "\\bfor\\b" << "\\bwhile\\b"
                    << "\\bif\\b" << "\\belse\\b" << "\\bvoid\\b" << "\\bunsigned\\b"
                    << "\\btrue\\b" << "\\bfalse\\b" << "\\breturn\\b" << "\\bstruct\\b"
                    << "\\bbool\\b" << "\\buchar\\b" << "\\bushort\\b" << "\\bint\\b" << "\\buint\\b" << "\\bulong\\b" << "\\bhalf\\b" << "\\bquad\\b" << "\\bfloat\\b"
                    << "\\bchar2\\b" << "\\bchar3\\b" << "\\bchar4\\b" << "\\bchar8\\b" << "\\bchar16\\b"
                    << "\\buchar2\\b" << "\\buchar3\\b" << "\\buchar4\\b" << "\\buchar8\\b" << "\\buchar16\\b"
                    << "\\bshort2\\b" << "\\bshort3\\b" << "\\bshort4\\b" << "\\bshort8\\b" << "\\bshort16\\b"
                    << "\\bushort2\\b" << "\\bushort3\\b" << "\\bushort4\\b" << "\\bushort8\\b" << "\\bushort16\\b"
                    << "\\bint2\\b" << "\\bint3\\b" << "\\bint4\\b" << "\\bint8\\b" << "\\bint16\\b"
                    << "\\buint2\\b" << "\\buint3\\b" << "\\buint4\\b" << "\\buint8\\b" << "\\buint16\\b"
                    << "\\blong2\\b" << "\\blong3\\b" << "\\blong4\\b" << "\\blong8\\b" << "\\blong16\\b"
                    << "\\bulong2\\b" << "\\bulong3\\b" << "\\bulong4\\b" << "\\bulong8\\b" << "\\bulong16\\b"
                    << "\\bfloat2\\b" << "\\bfloat3\\b" << "\\bfloat4\\b" << "\\bfloat8\\b" << "\\bfloat16\\b"
                    << "\\bdouble2\\b" << "\\bdouble3\\b" << "\\bdouble4\\b" << "\\bdouble8\\b" << "\\bdouble16\\b"
                    << "\\bhalf2\\b" << "\\bhalf3\\b" << "\\bhalf4\\b" << "\\bhalf8\\b" << "\\bhalf16\\b"
                    << "\\bptrdiff_t\\b" << "\\bintptr_t\\b" << "\\buintptr_t\\b"
                    << "\\bimage2d_t\\b" << "\\bimage3d_t\\b" << "\\bsampler_t\\b" << "\\bevent_t\\b"
                    << "\\bbool2\\b" << "\\bbool3\\b" << "\\bbool4\\b" << "\\bbool8\\b" << "\\bbool16\\b"
                    << "\\bquad2\\b" << "\\bquad3\\b" << "\\bquad4\\b" << "\\bquad8\\b" << "\\bquad16\\b"
                    << "\\b_cl_platform_id\\b" << "\\b_cl_device_id\\b" << "\\b_cl_context\\b" << "\\b_cl_command_queue\\b"
                    << "\\b_cl_mem\\b" << "\\b_cl_program\\b" << "\\b_cl_kernel\\b" << "\\b_cl_event\\b" << "\\b_cl_sampler\\b"
                    << "\\bcl_image_format\\b"
                    << "\\b(float|double)(2|3|4|8|16)x(2|3|4|8|16)\\b"
                    ;

    foreach(QString pattern, keywordPatterns) {
        rule.pattern = QRegExp(pattern);
        rule.format = m_keywordFormat;
        m_rules.append(rule);
    }

    // setup constants
    m_keywordFormat.setForeground(Qt::magenta);
    QStringList constantPatterns;

    constantPatterns << "\\bMAXFLOAT\\b" << "\\bHUGE_VALF\\b" << "\\bINFINITY\\b" << "\\bNAN\\b"
                     << "\\bFLT_DIG\\b" << "\\bFLT_MANT_DIG\\b" << "\\bFLT_MAX_10_EXP\\b"
                     << "\\bFLT_MAX_EXP\\b" << "\\bFLT_MIN_10_EXP\\b" << "\\bFLT_MIN_EXP\\b"
                     << "\\bFLT_RADIX\\b" << "\\bFLT_MAX\\b" << "\\bFLT_MIN\\b"
                     << "\\bFLT_EPSILON\\b" << "\\bCHAR_BIT\\b" << "\\bCHAR_MAX\\b"
                     << "\\bCHAR_MIN\\b" << "\\bINT_MIN\\b" << "\\bINT_MAX\\b" << "\\bLONG_MAX\\b"
                     << "\\bLONG_MIN\\b" << "\\bSCHAR_MAX\\b" << "\\bSCHAR_MIN\\b"
                     << "\\bSHRT_MAX\\b" << "\\bSHRT_MIN\\b" << "\\bUCHAR_MAX\\b"
                     << "\\bUCHAR_MIN\\b" << "\\bUINT_MAX\\b" << "\\bULONG_MAX\\b"
                     << "\\bDBL_DIG\\b" << "\\bDBL_MANT_DIG\\b" << "\\bDBL_MAX_10_EXP\\b"
                     << "\\bDBL_MIN_10_EXP\\b" << "\\bDBL_MIN_EXP\\b" << "\\bDBL_MAX\\b"
                     << "\\bDBL_MIN\\b" << "\\bDBL_EPSILON\\b" << "\\bM_E\\b" << "\\bM_LOG2E\\b"
                     << "\\bM_LOG10E\\b" << "\\bM_LN2\\b" << "\\bM_LN10\\b" << "\\bM_PI\\b"
                     << "\\bM_PI2\\b" << "\\bM_PI4\\b" << "\\bM_1_PI\\b" << "\\bM_2_PI\\b"
                     << "\\bM_2_SQRTPI\\b" << "\\bM_SQRT2\\b" << "\\bM_SQRT1_2\\b"
                     << "\\bCLK_NORMALIZED_COORDS_TRUE\\b" << "\\bCLK_NORMALIZED_COORDS_FALSE\\b"
                     << "\\bCLK_ADDRESS_REPEAT\\b" << "\\bCLK_ADDRESS_CLAMP_TO_EDGE\\b"
                     << "\\bCLK_ADDRESS_CLAMP\\b" << "\\bCL_INTENSITY\\b" << "\\bCL_RA\\b"
                     << "\\bCL_ARGB\\b" << "\\bCL_BGRA\\b" << "\\bCL_RGBA\\b" << "\\bCL_R\\b"
                     << "\\bCL_RG\\b" << "\\bCL_RGB\\b" << "\\bCL_RGx\\b" << "\\bCL_RGBx\\b"
                     << "\\bCL_Rx\\b" << "\\bCL_A\\b" << "\\bCL_LUMINANCE\\b"
                     << "\\bCL_SNORM_INT8\\b" << "\\bCL_SNORM_INT16\\b" << "\\bCL_UNORM_INT8\\b"
                     << "\\bCL_UNORM_INT16\\b" << "\\bCL_UNORM_SHORT_565\\b"
                     << "\\bCL_UNORM_SHORT_555\\b" << "\\bCL_UNORM_INT_101010\\b"
                     << "\\bCL_SIGNED_INT8\\b" << "\\bCL_SIGNED_INT16\\b"
                     << "\\bCL_SIGNED_INT32\\b" << "\\bCL_UNSIGNED_INT8\\b"
                     << "\\bCL_UNSIGNED_INT16\\b" << "\\bCL_UNSIGNED_INT32\\b"
                     << "\\bCL_HALF_FLOAT\\b" << "\\bCL_FLOAT\\b" << "\\bCLK_ADDRESS_NONE\\b"
                     << "\\bCLK_FILTER_NEAREST\\b" << "\\bCLK_FILTER_LINEAR\\b"
                     << "\\bCLK_GLOBAL_MEM_FENCE\\b" << "\\bCLK_LOCAL_MEM_FENCE\\b"
                     ;

    foreach(QString pattern, constantPatterns) {
        rule.pattern = QRegExp(pattern);
        rule.format = m_keywordFormat;
        m_rules.append(rule);
    }
}

void OpenCLHighlighter::setupReservedKeywords() {
    //highlightRule_t rule;
    //QStringList keywordPatterns;

    //m_reservedKeywordFormat.setForeground(Qt::red);
    //m_reservedKeywordFormat.setFontWeight(QFont::Bold);

    ////
    //// reserved keywords
    ////
    //keywordPatterns << "\\basm\\b" << "\\bclass\\b" << "\\bunion\\b"
                    //<< "\\benum\\b" << "\\btypedef\\b" << "\\btemplate\\b"
                    //<< "\\bthis\\b" << "\\bpacked\\b" << "\\bgoto\\b"
                    //<< "\\bswitch\\b" << "\\bdefault\\b" << "\\binline\\b"
                    //<< "\\bnoinline\\b" << "\\bvolatile\\b" << "\\bpublic\\b"
                    //<< "\\bstatic\\b" << "\\bextern\\n" << "\\bexternal\\b"
                    //<< "\\binterface\\b" << "\\blong\\b" << "\\bshort\\b"
                    //<< "\\bdouble\\b" << "\\bhalf\\b" << "\\bfixed\\b"
                    //<< "\\bunsigned\\b" << "\\blowp\\b" << "\\bmediump\\b"
                    //<< "\\bhighp\\b" << "\\bprecision\\b"
                    //<< "\\binput\\b" << "\\boutput\\b"
                    //<< "\\bhvec2\\b" << "\\bhvec3\\b" << "\\bhvec4\\b"
                    //<< "\\bdvec2\\b" << "\\bdvec3\\b" << "\\bdvec4\\b"
                    //<< "\\bfvec2\\b" << "\\bfvec3\\b" << "\\bfvec4\\b"
                    //<< "\\bsampler2DRect\\b" << "\\bsampler3DRect\\b" << "\\bsampler2DRectShadow\\b"
                    //<< "\\bsizeof\\b" << "\\bcast\\b" << "\\bnamespace\\b" << "\\busing\\b";

    //foreach(QString pattern, keywordPatterns) {
        //rule.pattern = QRegExp(pattern);
        //rule.format = m_reservedKeywordFormat;
        //m_rules.append(rule);
    //}
}

} // namespace
