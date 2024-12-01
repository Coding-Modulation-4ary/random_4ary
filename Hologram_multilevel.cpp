#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "uniform&gaussianRNG.h"
#include "chi_46proposal_modulation_code.h"

/*******************	Definition of Constants	**************************/
const double PI = 3.14159;						////  Value of pi
const double DELTA = 0.000005;						//// Defines the range of variation for the coefficients of the Equalizer. // Độ nhạy của bộ cân bằng tín hiệu (Equalizer). Giá trị nhỏ giúp tăng độ chính xác.

//const double DELTA				=	0.000005;						//// Defines the range of variation for the coefficients of the Equalizer
																	//// If the value is large, the variation of the coefficients will be more stable,
								//// but it may fail to find the optimal value. If the variation is small, it may converge more precisely.
const int Training_Block = 3; // Số khối dữ liệu dùng để huấn luyện mô hình.
const int M_ary = 4;								//// the number of total symbols (2 ary = 1but (binary), 4 ary = 2bit, 3 ary = 1.5bit)
const int BLOCK = 5;      						//// the total block. Tổng số khối dữ liệu để xử lý.
const int Page_Size = 1024;							//// Block size of the hologram. Kích thước của mỗi khối dữ liệu.
const int Data_Size = Page_Size;					//// Kích thước dữ liệu tương tự như Page_Size.
const double SigmaB = 1.0;							////  Noise level of the hologram. Mức nhiễu (Noise level) trong tín hiệu holographic.
const int SELL = 5;								//// pixel center==0, Sell==odd.  Pixel trung tâm (giá trị lẻ).
const int Unit_Sell = 11;								//// The one block of one sell divided. Một khối dữ liệu được chia thành 11 phần nhỏ.

const int Equalizer_Length = 5;								//// must set odd. Độ dài của bộ cân bằng tín hiệu (phải là số lẻ để cân bằng hai phía).
const int Pridictor_Length = 3;								////  Length of the predictor. Độ dài của bộ dự đoán tín hiệu.
const int Adaptive_Length = 5;								//// Adaptive length. Độ dài của bộ điều chỉnh tín hiệu thích nghi (Adaptive Equalizer).

const double Mis_const = 1;								//// Specifies the constant for Mis-alignment. Chỉ định giá trị hằng cho độ lệch (Mis-alignment)
//// 0: Random, 1: const value (Mis_alignment_x, Mis_alignment_y)
const double Mis_alignment_dB = 30;  // Mức nhiễu (dB) do sai lệch.
const double Mis_alignment_x = 0.0;  // Độ lệch trên trục X.
const double Mis_alignment_y = 0.0;  // Độ lệch trên trục Y.
const int Target_Order = 3;								//// Order of the PR Target filter. Must be an odd number. // Bậc của bộ lọc mục tiêu PR (phải là số lẻ).
const int depth = 29;								//// viterbi constraint // Độ sâu ràng buộc của thuật toán Viterbi.

const int noise_free = 1;								//// 0 : noise free, 1: noise exist
const double START_SNR = 10;//int((SigmaB+0.0001)/0.05-36);  // Giá trị SNR ban đầu.
const double END_SNR = START_SNR; // Giá trị SNR kết thúc.
const double SNR_INTERVAL = 1; // Bước nhảy của SNR giữa các lần lặp.

const char* mis[3] = { "mis0", "mis10", "mis20" }; // Các mức sai lệch khác nhau để kiểm tra.
const int checkHist = 0; //1, lay hist // 1: Lưu histogram, 0: Không lưu.
const int Encoder = 11; // Phương pháp mã hóa: 
//// 0: random, 10: 24mc_chi, 11: 46mc_chi
const char* s_encoder[10] = { "ran","prof_park","kim_3_4","sun_park_69","kim_2_3","seungmin_69","gg_grad_23_modulation_kim","gg_grad_23_modulation_prof","gg_dmin2","gg_dmin2_state" };  // Các phương pháp mã hóa tín hiệu.
const int DRAWING_BER_CURVE = 1;				//// 0: float, 1: int

double POWER = 1;							//// To maintain a consistent power level. // Duy trì mức công suất tín hiệu nhất quán.

const int DETECTOR = 0;
//// 0: 1D detector, // Bộ phát hiện tín hiệu. 0: 1D detector.

const int Modulation_Decoder = 0;	//// For ECC 4/6 codes. 0: Viterbi, 1: SOVA // Giải mã tín hiệu. 0: Viterbi, 1: SOVA.
const int Iteration_Det_Dem = 1;	//// # of the iteration between detecter(2D SOVA) and demodulation (SOVA for demodulation)  // Số lần lặp giữa phát hiện tín hiệu và giải mã.
const char* s_detecter[8] = { "1D_Viterbi","2D_NPSOVA","2D NPSOVA Hard","Only AWGN","Adaptive 2D NPSOVA","proposed iterative 2D SOVA","2D Log-MAP","gg_grad_A" };  // Các phương pháp phát hiện tín hiệu.
const int Quantization = 0;								//// 1,2,3,4 : Quantization, else: None Quantization, // 1,2,3,4: Có lượng tử hóa. Khác: Không lượng tử hóa.

const int NoiseFilter = 0;								//// 0: No noisefilter,  1: use noisefilter  // 0: Không dùng bộ lọc nhiễu, 1: Dùng bộ lọc nhiễu.
const int Adaptive_yesornot = 0;								//// 0: No adaptive,  1: use adaptive  // 0: Không dùng Adaptive Equalizer, 1: Dùng Adaptive Equalizer.
const int AdaptiveNF_yesornot = 0;								//// 0: No adaptive Nosefilter,  1: use adaptive Nosefilter // 0: Không dùng bộ lọc nhiễu thích nghi, 1: Có dùng.

const int M_number = 4;								//// 0: no adapt M-alorithm, 1~ : adapt M-alorithm number  // Số giá trị cho thuật toán M-algorithm.
const int Constant_period = 0;								//// In M-algirhtm, not use M-algorithm from Constant_period to 0 reversly // Không sử dụng M-algorithm từ Constant_period đến 0 ngược lại.
const int lee_proposed_method = 0;	//// 0: 1st-8state, 1: 2nd-8state, 2: 3rd-0state, 3: 4th-16state  // Phương pháp đề xuất: 
const int lee_proposed_constraint = 30;	//// 30, 100, 500, 2000 // Ràng buộc: 30, 100, 500, 2000.

const int SCALE_DEMODUL = 10;								//// During de-modulating, the number applied Viterbi decoder  // Số lần áp dụng giải mã trong quá trình de-modulation.

long k_noise = 1;   // Hệ số nhiễu trong quá trình xử lý tín hiệu.

//////////////////////// In Viterbi, cases with 03 and 30 states are excluded ///////////////////////////////////////
const int gg_pr_viterbi = 0;		 // 1: Phương pháp Viterbi đề xuất, khác: Viterbi thông thường.					///     1 : proposed, else : conventional
const int gg_reduce_last = 0;   // Tùy chọn giảm số trạng thái cuối.
int ggg = 0;  // Biến tạm thời (debug hoặc điều khiển trạng thái).
const int modul_method = 1;							///		0 : state 8, 1 : state 16  // Phương pháp điều chế: 0: state 8, 1: state 16.

/*
 * Performs 2D quantization on the input data matrix.
 * Parameters:
 * - data: Input 2D matrix of data.
 * - outdata: Output 2D matrix after quantization.
 * - Page_Size: Number of elements in each page.
 * - max: Maximum value in the data range.
 * - min: Minimum value in the data range.
 * - quantization_number: Number of levels for quantization.
 */
void TwoD_Quantization(double** data, double** outdata, int Page_Size, double max, double min, int quantization_number);
/*
 * Simulates the holographic channel by applying noise and distortion.
 * Parameters:
 * - h: Output 2D matrix representing the channel response.
 * - SigmaB: Noise level for the channel.
 * - Mis_alignment_sigma: Misalignment noise level.
 * - SELL: Size of pixel center (odd number).
 * - Unit_Sell: Number of subdivisions in one sell.
 */
void HoloChannel(double** h, double SigmaB, double Mis_alignment_sigma, int SELL, int Unit_Sell);
/*
 * Sets the coefficients for the PR Target filter.
 * Parameters:
 * - Target_Coef: Output 2D array of filter coefficients.
 * - target_order: Order of the filter (must be odd).
 */
void PR_Target_set(double** Target_Coef, int target_order);
/*
 * Performs convolution with PRML (Partial Response Maximum Likelihood).
 * Parameters:
 * - trellis: 2D array representing the trellis structure.
 * - Target_Coef: Coefficients of the PR Target filter.
 * - Target_Order: Order of the PR Target filter.
 */
void PRML_Convolution(double** trellis, double* Target_Coef, int Target_Order);
/*
 * Simulates the output of the holographic channel.
 * Parameters:
 * - outCH: Output 2D array after channel processing.
 * - h: Channel response matrix.
 * - inCH: Input 2D array representing the channel input.
 * - Page_Size: Size of the page being processed.
 * - SELL: Size of pixel center (odd number).
 * - sigma: Noise level in the channel.
 */
void Hologram_ChannelOut(double** outCH, double** h, int** inCH, int Page_Size, int SELL, double sigma);
/*
 * Applies an equalizer in the X direction for holographic data.
 * Parameters:
 * - outEQ: Output equalized data.
 * - eq_coef: Equalizer coefficients.
 * - inEQ: Input data to the equalizer.
 * - Page_Size: Size of the data page.
 * - SELL: Pixel center size.
 * - Equalizer_Length: Length of the equalizer filter.
 * - Training: Whether training mode is enabled.
 * - inCH: Input channel data.
 * - Target_Coef: Coefficients of the PR Target filter.
 * - Target_Order: Order of the PR Target filter.
 */
void Hologram_EqualizerOut_x(double** outEQ, double** eq_coef, double** inEQ, int Page_Size, int SELL, int Equalizer_Length, int Training, int** inCH, double** Target_Coef, int Target_Order);
/*
 * Applies an equalizer in the Y direction for holographic data.
 * (Parameters are similar to Hologram_EqualizerOut_x.)
 */
void Hologram_EqualizerOut_y(double** outEQ, double** eq_coef, double** inEQ, int Page_Size, int SELL, int Equalizer_Length, int Training, int** inCH, double** Target_Coef, int Target_Order);
/*
 * Generates predictor coefficients in the X direction for holographic data.
 * Parameters:
 * - Pridictor_Coef_x: Output predictor coefficients.
 * - inEQ: Input equalized data.
 * - Page_Size: Size of the data page.
 * - SELL: Pixel center size.
 * - Equalizer_Length: Length of the equalizer filter.
 * - inCH: Input channel data.
 * - Target_Coef: Coefficients of the PR Target filter.
 * - Target_Order: Order of the PR Target filter.
 * - Pridictor_Length: Length of the predictor filter.
 */

void NPML_Pridictor_x(double* Pridictor_Coef_x, double** inEQ, int Page_Size, int SELL, int Equalizer_Length, int** inCH, double** Target_Coef, int Target_Order, int Pridictor_Length);
/*
 * Similar to NPML_Pridictor_x but operates in the Y direction.
 */
void NPML_Pridictor_y(double* Pridictor_Coef_y, double** inEQ, int Page_Size, int SELL, int Equalizer_Length, int** inCH, double** Target_Coef, int Target_Order, int Pridictor_Length);
/*
 * Performs SOVA decoding for holographic data in the X direction.
 * Parameters:
 * - out_SOVA: Output data after SOVA decoding.
 * - eq_out: Equalized input data.
 * - EI: Extrinsic information for decoding.
 * - trellis: Trellis structure for decoding.
 * - Pridictor_Coef: Predictor coefficients.
 * - Pridictor_Length: Length of the predictor.
 * - Target_Order: Order of the PR Target filter.
 * - Page_Size: Size of the data page.
 * - depth: Depth constraint of the trellis.
 */
void Hologram_SOVA_x(double** out_SOVA, double** eq_out, double** EI, double** trellis, double* Pridictor_Coef, int Pridictor_Length, int Target_Order, int Page_Size, int depth);
/*
 * Similar to Hologram_SOVA_x but operates in the Y direction.
 */
void Hologram_SOVA_y(double** out_SOVA, double** eq_out, double** EI, double** trellis, double* Pridictor_Coef, int Pridictor_Length, int Target_Order, int Page_Size, int depth);
/*
 * Generates interleaver and de-interleaver tables for holographic data.
 * Parameters:
 * - Int: Interleaver table.
 * - DeInt: De-interleaver table.
 * - nSize: Size of the table.
 */
void MakeInterleaverFile(int* Int, int* DeInt, int nSize);
/*
 * Applies interleaving to holographic data (integer version).
 * Parameters:
 * - Page: Input page data.
 * - Inter_Page: Output interleaved page data.
 * - Inter_table: Interleaver table.
 * - Page_Size: Size of the page.
 */
void Interleaver_for_Holo(int** Page, int** Inter_Page, int* Inter_table, int Page_Size);
/*
 * Applies interleaving to holographic data (double version).
 */
void Interleaver_for_Holo(double** Page, double** Inter_Page, int* Inter_table, int Page_Size);
/*
 * Performs PRML convolution for holographic data with GG Gradient 23.
 */
void PRML_Convolution_gg_grad_23(double** trellis, double* Target_Coef, int Target_Order, int input_size, int* input);
/*
 * Similar to PRML_Convolution_gg_grad_23 but applies to specific modular trellises.
 */
void PRML_Convolution_for_modul_gg_grad_23(int*** trellis_modul);
void PRML_Convolution_for_modul_gg_grad_23_prof(int*** trellis_modul);
void PRML_Convolution_for_modul_gg_grad_23_prof_method0(int**** trellis_modul);
void PRML_Convolution_for_modul_gg_dmin2(int*** trellis_modul);
void PRML_Convolution_for_modul_gg_dmin2_method0(int**** trellis_modul);

/*
 * Performs SOVA decoding with additional trellis data for holographic data in the X direction.
 * Parameters are similar to Hologram_SOVA_x but include an additional trellisA parameter.
 */
void Hologram_SOVA_x_for_gg_grad(double** out_SOVA, double** eq_out, double** EI, double** trellis, double** trellisA, double* Pridictor_Coef, int Pridictor_Length, int Target_Order, int Page_Size, int depth);


void main()
{
	clock_t c_time;
	c_time = clock(); // c_time Check Start 

	int i, j, k; // Biến vòng lặp

	int A_input[2] = { 0,1 }; // 2 biến ban đầu là 0 và 1
	int BC_input[4] = { 0,1,2,3 }; // 4-ary 
	int total_states_method0 = 0; // Đếm tổng số trạng thái của hệ thống hoặc thuật toán được cấu hình ở method0.

	double** trellis_A; // Biến này sẽ được sử dụng để lưu cấu trúc trellis cho PRML hoặc SOVA decoding. Cấu trúc trellis là một phần quan trọng trong giải mã tín hiệu modulation code. Dữ liệu này có thể chứa các giá trị trạng thái hoặc xác suất chuyển tiếp giữa các trạng thái.


	int file_see = 0; // Biến này có thể được sử dụng để đánh dấu trạng thái của tệp hoặc để kiểm tra xem tệp có được xử lý (hoặc mở) thành công không.

	FILE* OUT, * BRIEF_RESULT, * BRIEF_RESULT_MIS, * fber, * CHOUT[M_ary], * gg_fp_en, * gg_fp_de, * gg_fp_pattern; // Các biến này là con trỏ tới tệp được sử dụng để ghi dữ liệu kết quả hoặc dữ liệu trung gian trong quá trình thực thi chương trình.

	fber = fopen("test.dat", "w"); // Tệp dùng để ghi dữ liệu thử nghiệm. Dữ liệu trong tệp này có thể là các chỉ số hiệu suất hoặc tín hiệu kiểm tra.
	BRIEF_RESULT = fopen("breif_result.m", "a"); // Lưu các kết quả chung.
	BRIEF_RESULT_MIS = fopen("breif_result_mis.m", "a"); // Lưu các kết quả liên quan đến sai lệch (misalignment).
	// CHOUT[M_ary]:  Một mảng con trỏ tệp, số lượng tệp phụ thuộc vào giá trị M_ary (4 trong trường hợp 4-ary modulation). Lưu tín hiệu đầu ra kênh ứng với các trạng thái modulation khác nhau.
	gg_fp_en = fopen("encoded.dat", "a");
	gg_fp_de = fopen("decoded.dat", "a");
	gg_fp_pattern = fopen("pattern.dat", "a");
	// Các mảng ký tự này được dùng để lưu đường dẫn hoặc tên tệp liên quan đến các mức sai lệch (misalignment).
	// Kích thước 255 đảm bảo đủ lớn để chứa tên tệp dài.
	char nameMis0[255];
	char nameMis1[255];
	char nameMis2[255];
	char nameMis3[255];
	char nameMis4[255];
	char nameMis5[255];
	// Các con trỏ tệp này lưu dữ liệu đầu ra từ các bước khác nhau của quá trình xử lý tín hiệu.
	// Chúng được tổ chức theo loại dữ liệu:
	FILE* fber_outCH, * fber_outCH_Mary; // Lưu dữ liệu tín hiệu từ kênh (Channel Output).
	FILE* fber_outSOVA, * fber_outSOVA_Mary; // Lưu dữ liệu sau bước giải mã SOVA.
	FILE* fber_outDEM, * fber_outDEM_Mary; //  Lưu dữ liệu từ bộ giải điều chế (Demodulator).

	// Đoạn mã này tạo ra tên tệp động dựa trên giá trị của Mis_alignment_x, sau đó mở các tệp tương ứng để ghi dữ liệu đầu ra.
	// Mis_alignment_x đại diện cho mức độ sai lệch trong trục X và được sử dụng để phân loại dữ liệu đầu ra theo các trường hợp sai lệch khác nhau.
	if (Mis_alignment_x == 0)
	{
		sprintf(nameMis0, "%s_fber_outCH.txt", mis[0]);
		sprintf(nameMis1, "%s_fber_outCH_Mary.txt", mis[0]);
		sprintf(nameMis2, "%s_fber_outSOVA.txt", mis[0]);
		sprintf(nameMis3, "%s_fber_outSOVA_Mary.txt", mis[0]);
		sprintf(nameMis4, "%s_fber_outDEM.txt", mis[0]);
		sprintf(nameMis5, "%s_fber_outDEM_Mary.txt", mis[0]);

		fber_outCH = fopen(nameMis0, "a");
		fber_outCH_Mary = fopen(nameMis1, "a");

		fber_outSOVA = fopen(nameMis2, "a");
		fber_outSOVA_Mary = fopen(nameMis3, "a");

		fber_outDEM = fopen(nameMis4, "a");
		fber_outDEM_Mary = fopen(nameMis5, "a");

	}
	else if (Mis_alignment_x == 0.1)
	{
		sprintf(nameMis0, "%s_fber_outCH.txt", mis[1]);
		sprintf(nameMis1, "%s_fber_outCH_Mary.txt", mis[1]);
		sprintf(nameMis2, "%s_fber_outSOVA.txt", mis[1]);
		sprintf(nameMis3, "%s_fber_outSOVA_Mary.txt", mis[1]);
		sprintf(nameMis4, "%s_fber_outDEM.txt", mis[1]);
		sprintf(nameMis5, "%s_fber_outDEM_Mary.txt", mis[1]);

		fber_outCH = fopen(nameMis0, "a");
		fber_outCH_Mary = fopen(nameMis1, "a");

		fber_outSOVA = fopen(nameMis2, "a");
		fber_outSOVA_Mary = fopen(nameMis3, "a");

		fber_outDEM = fopen(nameMis4, "a");
		fber_outDEM_Mary = fopen(nameMis5, "a");
	}
	else
	{
		sprintf(nameMis0, "%s_fber_outCH.txt", mis[2]);
		sprintf(nameMis1, "%s_fber_outCH_Mary.txt", mis[2]);
		sprintf(nameMis2, "%s_fber_outSOVA.txt", mis[2]);
		sprintf(nameMis3, "%s_fber_outSOVA_Mary.txt", mis[2]);
		sprintf(nameMis4, "%s_fber_outDEM.txt", mis[2]);
		sprintf(nameMis5, "%s_fber_outDEM_Mary.txt", mis[2]);

		fber_outCH = fopen(nameMis0, "a");
		fber_outCH_Mary = fopen(nameMis1, "a");

		fber_outSOVA = fopen(nameMis2, "a");
		fber_outSOVA_Mary = fopen(nameMis3, "a");

		fber_outDEM = fopen(nameMis4, "a");
		fber_outDEM_Mary = fopen(nameMis5, "a");
	}


	unsigned long count = 0; // Biến đếm tổng số lần thực thi hoặc số lượng tín hiệu được xử lý trong toàn bộ chương trình.
	unsigned long block_count = 0; // Biến đếm số khối dữ liệu (blocks) đã được xử lý.
	unsigned long count_M_ary = 0; //  Biến đếm số tín hiệu được xử lý trong modulation code dạng M-ary (ví dụ: 4-ary hoặc 16-ary modulation).
	unsigned long count_modul = 0;// Đếm số tín hiệu đã được mã hóa (modulated).
	unsigned long count_modul_M_ary = 0; // Đếm số tín hiệu đã được mã hóa và phân loại theo modulation dạng M-ary.
	unsigned long error = 0; //  Biến đếm tổng số lỗi phát sinh trong quá trình xử lý tín hiệu.
	unsigned long error_M_ary = 0; //  Đếm số lỗi phát sinh trong các tín hiệu được mã hóa bằng modulation dạng M-ary.
	unsigned long error_eq = 0; // Đếm số lỗi phát sinh sau bước cân bằng tín hiệu (Equalizer).
	unsigned long error_eq_other = 0; // Đếm lỗi phát sinh trong các trường hợp khác không nằm trong bước cân bằng tín hiệu.
	unsigned long error_viterbi[Iteration_Det_Dem * 2] = { 0, }; // Lưu lỗi phát sinh trong thuật toán giải mã Viterbi. Kích thước của mảng phụ thuộc vào Iteration_Det_Dem, biểu thị số lần lặp giữa phát hiện tín hiệu và giải mã.
	unsigned long error_viterbi_M_ary[Iteration_Det_Dem * 2] = { 0, }; // Lưu lỗi giải mã Viterbi cho tín hiệu modulation dạng M-ary.
	unsigned long error_viterbi_M_ary_gg[6] = { 0, }; // Mảng lưu lỗi giải mã Viterbi liên quan đến các thuật toán hoặc cấu hình modulation cụ thể (được biểu thị bằng "gg").
	unsigned long error_demodul[Iteration_Det_Dem * 2] = { 0, }; // Lưu lỗi phát sinh trong bước giải điều chế.
	unsigned long error_demodul_M_ary[Iteration_Det_Dem * 2] = { 0, }; // Lưu lỗi giải điều chế cho tín hiệu modulation dạng M-ary.
	unsigned long error_demodul_check[Iteration_Det_Dem] = { 0, }; // Lưu lỗi giải điều chế cho tín hiệu modulation dạng M-ary.

	int temp_count[M_ary] = { 0, }; // Mảng lưu số lượng tín hiệu thuộc từng trạng thái trong modulation dạng M-ary. Ví dụ: Nếu M_ary = 4, mảng này sẽ lưu số lượng tín hiệu cho các trạng thái {0, 1, 2, 3}.
	double temp_sum[M_ary] = { 0, }; // Mảng lưu tổng giá trị tín hiệu tương ứng với từng trạng thái. Dùng để tính trung bình hoặc phân tích các trạng thái của modulation code. Dùng để tính trung bình hoặc phân tích các trạng thái của modulation code.

	double mis_sigma = 0.0; // Lưu giá trị sai lệch (misalignment noise) được tính toán từ SNR hoặc các thông số khác.
	double sigma = 0.0;										//// Using SNR_DB to calculate and store the value of sigma. Độ lệch chuẩn (standard deviation) của nhiễu Gaussian trong kênh truyền. Giá trị này được tính toán từ SNR (Signal-to-Noise Ratio).
	double snr;												//// Using SNR_DB to calculate and store the value of SNR. Tỉ số tín hiệu trên nhiễu (SNR), được lưu dưới dạng số thực (linear SNR).
	double eqout_total_x = 0, eqout_average_x = 0, eqout_total_y = 0, eqout_average_y = 0; // Tổng và trung bình của tín hiệu đầu ra sau khi qua bộ cân bằng (Equalizer) theo trục X. Tổng và trung bình của tín hiệu đầu ra sau khi qua bộ cân bằng (Equalizer) theo trục Y.
	double chout_average[M_ary - 1] = { 0, }; // Mảng lưu giá trị trung bình của tín hiệu đầu ra từ kênh truyền (channel output) cho các trạng thái modulation code (trừ trạng thái cuối cùng).
	char strH[255]; // Mảng ký tự để lưu tên hoặc đường dẫn của tệp liên quan đến kênh truyền (channel).
	// Kích thước 255 đảm bảo đủ lớn để chứa các chuỗi dài.

	double** h;												//// Storing the channel. Ma trận lưu thông tin của kênh truyền (channel response). Được sử dụng để mô phỏng quá trình tín hiệu đi qua kênh và thêm nhiễu.
	double** Target_Coef;									//// Setting the coefficients of the PR Target. Ma trận lưu hệ số của bộ lọc mục tiêu PR (Partial Response Target). Hệ số này được sử dụng trong các thuật toán như PRML Convolution hoặc Viterbi.

	double** trellis; // Ma trận lưu cấu trúc trellis dùng cho giải mã. Trellis đại diện cho các trạng thái và xác suất chuyển tiếp trong thuật toán SOVA hoặc Viterbi.

	// Allocate memory for the Target_Coef matrix (Target_Order x Target_Order)
	// This matrix stores the coefficients for the PR Target filter.
	// Sau khi được tạo, nó được truyền vào các hàm xử lý như convolution, equalization, hoặc decoding, nơi các giá trị trong ma trận sẽ được sử dụng hoặc cập nhật.
	Target_Coef = (double**)malloc(Target_Order * sizeof(double*)); // Allocate rows
	for (i = 0; i < Target_Order; i++) {
		Target_Coef[i] = (double*)malloc(Target_Order * sizeof(double)); // Allocate columns for each row
	}
	// Link docs cụ thể về vai trò của Target_Coef trong bộ lọc PR Target: https://docs.google.com/document/d/1L9zRDMRuHSnI_ra8Q8xmUUons77RcviouuEtCtuT90Q/edit?tab=t.0

	////// For Viterbi and BCJR ///////
    ////// Allocate memory for trellis structure used in Viterbi and BCJR algorithms ///////
	trellis = (double**)malloc((int)pow((double)M_ary, Target_Order - 1) * sizeof(double*));
	for (i = 0; i < pow((double)M_ary, Target_Order - 1); i++)
		trellis[i] = (double*)calloc(M_ary * 2, sizeof(double)); // // Allocate columns for each row and initialize to 0

	////// For Viterbi and BCJR : gukhui's graduation ///////
	trellis_A = (double**)malloc((int)pow((double)M_ary - 2, Target_Order - 1) * sizeof(double*));
	for (i = 0; i < pow((double)M_ary - 2, Target_Order - 1); i++)
		trellis_A[i] = (double*)calloc((M_ary - 2) * 2, sizeof(double));

	int*** trellis_modul = (int***)malloc((int)pow((double)M_ary, 2) * sizeof(int**));
	for (i = 0; i < (int)pow((double)M_ary, 2); i++)
		trellis_modul[i] = (int**)malloc((int)pow((double)M_ary, 2) * sizeof(int*));
	for (i = 0; i < (int)pow((double)M_ary, 2); i++)
		for (j = 0; j < (int)pow((double)M_ary, 2); j++)
			trellis_modul[i][j] = (int*)calloc(3, sizeof(int));

	int*** trellis_modul_prof = (int***)malloc((int)pow((double)M_ary, 2) * sizeof(int**));
	for (i = 0; i < (int)pow((double)M_ary, 2); i++)
		trellis_modul_prof[i] = (int**)malloc((int)pow((double)M_ary, 2) * sizeof(int*));
	for (i = 0; i < (int)pow((double)M_ary, 2); i++)
		for (j = 0; j < (int)pow((double)M_ary, 2); j++)
			trellis_modul_prof[i][j] = (int*)calloc(3, sizeof(int));

	total_states_method0 = ((int)pow((double)M_ary, 2)) / 2;
	int**** trellis_modul_prof_method0 = (int****)malloc(total_states_method0 * sizeof(int***));
	for (i = 0; i < total_states_method0; i++)
		trellis_modul_prof_method0[i] = (int***)malloc(total_states_method0 * sizeof(int**));
	for (i = 0; i < total_states_method0; i++)
		for (j = 0; j < total_states_method0; j++)
			trellis_modul_prof_method0[i][j] = (int**)malloc(2 * sizeof(int*));
	for (i = 0; i < total_states_method0; i++)
		for (j = 0; j < total_states_method0; j++)
			for (k = 0; k < 2; k++)
				trellis_modul_prof_method0[i][j][k] = (int*)calloc(3, sizeof(int));


	///// For channel modeling ///////
	h = (double**)malloc(SELL * sizeof(double*));
	for (i = 0; i < SELL; i++)
		h[i] = (double*)calloc(SELL, sizeof(double));

	///// For PR equalizer //////
	double** Equalizer_Coef_x = (double**)malloc(Equalizer_Length * sizeof(double*));	//// Storing the coefficients of the Equalizer.
	double** Equalizer_Coef_y = (double**)malloc(Equalizer_Length * sizeof(double*));	//// Storing the coefficients of the Equalizer.
	for (i = 0; i < Equalizer_Length; i++)
		Equalizer_Coef_x[i] = (double*)calloc(Equalizer_Length, sizeof(double)), Equalizer_Coef_y[i] = (double*)calloc(Equalizer_Length, sizeof(double));

	double* Pridictor_Coef_x = (double*)calloc(Pridictor_Length + 1, sizeof(double));	//// Storing the coefficients of the Equalizer.
	double* Pridictor_Coef_y = (double*)calloc(Pridictor_Length + 1, sizeof(double));	//// Storing the coefficients of the Equalizer.

	///// For modulation codes /////
	// For Park 2/3 4-ary modulation
	int* input_park_2_3_4ary = (int*)calloc(Data_Size * (Data_Size / 3) * 2, sizeof(int));  // Input data for Park 2/3 modulation
	int* output_park_2_3_4ary = (int*)calloc(Data_Size * (Data_Size / 3) * 2, sizeof(int)); // Output data for Park 2/3 modulation

	// For Kim 3/4 4-ary modulation
	int* input_kim_3_4_4ary = (int*)calloc((Data_Size / 2) * (Data_Size / 2) * 3, sizeof(int));  // Input data for Kim 3/4 modulation
	int* output_kim_3_4_4ary = (int*)calloc((Data_Size / 2) * (Data_Size / 2) * 3, sizeof(int)); // Output data for Kim 3/4 modulation

	// For Sun-Park 6/9 4-ary modulation
	int* input_sun_park_6_9_4ary = (int*)calloc((Data_Size / 3) * (Data_Size / 3) * 6, sizeof(int));  // Input data for Sun-Park 6/9 modulation
	int* output_sun_park_6_9_4ary = (int*)calloc((Data_Size / 3) * (Data_Size / 3) * 6, sizeof(int)); // Output data for Sun-Park 6/9 modulation

	// For Chi 4/6 4-ary modulation
	int* input_chi_46_4ary = (int*)calloc((Data_Size / 3) * (Data_Size / 2) * 4, sizeof(int));  // Input data for Chi 46 modulation
	int* output_chi_46_4ary = (int*)calloc((Data_Size / 3) * (Data_Size / 2) * 4, sizeof(int)); // Output data for Chi 46 modulation

	// For Chi 24 4-ary modulation
	int* input_chi_24_4ary = (int*)calloc(Data_Size * (Data_Size / 4) * 2, sizeof(int));  // Input data for Chi 24 modulation
	int* output_chi_24_4ary = (int*)calloc(Data_Size * (Data_Size / 4) * 2, sizeof(int)); // Output data for Chi 24 modulation


	///// Allocate memory for intermediate data and channel operations /////

// Memory for extrinsic information (EI) and LDPC output
	double** EI = (double**)malloc(Page_Size * sizeof(double*)); // Extrinsic Information matrix
	int** outLDPC = (int**)malloc(Page_Size * sizeof(int*));     // LDPC decoded output
	for (i = 0; i < Page_Size; i++) {
		outLDPC[i] = (int*)calloc(Page_Size, sizeof(int));       // Initialize LDPC output
		EI[i] = (double*)calloc(Page_Size, sizeof(double));      // Initialize extrinsic information
	}

	// Memory for LDPC binary input tín hiệu đầu vào
	int** inLDPC = (int**)malloc(Data_Size * sizeof(int*));      // Binary input for LDPC
	for (i = 0; i < Data_Size; i++) {
		inLDPC[i] = (int*)calloc(Data_Size * M_ary, sizeof(int)); // Allocate binary input space
	}

	// Memory for LDPC M-ary symbol input Lưu trữ tín hiệu đầu vào ở dạng ký hiệu M-ary.
	int** inLDPC_symbol = (int**)malloc(Data_Size * sizeof(int*)); // M-ary input for LDPC
	for (i = 0; i < Data_Size; i++) {
		inLDPC_symbol[i] = (int*)calloc(Data_Size, sizeof(int));   // Allocate symbol input space
	}

	// Memory for channel input Tín hiệu đầu vào cho kênh.
	int** inCH = (int**)malloc((Page_Size + (Equalizer_Length / 2) * 2) * sizeof(int*)); // Input to the channel
	for (i = 0; i < Page_Size + (Equalizer_Length / 2) * 2; i++) {
		inCH[i] = (int*)calloc(Page_Size + (Equalizer_Length / 2) * 2, sizeof(int));     // Initialize channel input
	}

	// Memory for channel output Tín hiệu đầu ra sau khi đi qua kênh.
	double** outCH = (double**)malloc((Page_Size + (Equalizer_Length / 2) * 2) * sizeof(double*)); // Output after passing through the channel

	//// Allocate memory for channel output matrix (outCH) ////
	for (i = 0; i < Page_Size + (Equalizer_Length / 2) * 2; i++)
		outCH[i] = (double*)calloc((Page_Size + (Equalizer_Length / 2) * 2), sizeof(double)); // Allocate and initialize channel output

	//// Allocate memory for equalizer outputs (X and Y directions) 
	double** outEQ_x = (double**)malloc(Page_Size * sizeof(double*));  // Equalizer output (X-direction)
	double** outEQ_y = (double**)malloc(Page_Size * sizeof(double*));  // Equalizer output (Y-direction)

	//// Allocate memory for Viterbi outputs ////
	int** n_outViterbi = (int**)malloc(Page_Size * sizeof(int*)); // Integer output from Viterbi
	double** d_outViterbi = (double**)malloc(Page_Size * sizeof(double)); // Double output from Viterbi
	double** d_outViterbi_x = (double**)malloc(Page_Size * sizeof(double));  // Double output (X-direction)
	double** d_outViterbi_y = (double**)malloc(Page_Size * sizeof(double));  // Double output (Y-direction)

	for (i = 0; i < Page_Size; i++) {
		outEQ_x[i] = (double*)calloc(Page_Size, sizeof(double));  // Allocate and initialize equalizer output (X-direction)
		outEQ_y[i] = (double*)calloc(Page_Size, sizeof(double));  // Allocate and initialize equalizer output (Y-direction)
		n_outViterbi[i] = (int*)calloc(Page_Size, sizeof(int));   // Allocate and initialize integer Viterbi output
		d_outViterbi[i] = (double*)calloc(Page_Size, sizeof(double)); // Allocate and initialize double Viterbi output
		d_outViterbi_x[i] = (double*)calloc(Page_Size, sizeof(double)); // Allocate and initialize double Viterbi output (X-direction)
		d_outViterbi_y[i] = (double*)calloc(Page_Size, sizeof(double)); // Allocate and initialize double Viterbi output (Y-direction)
	}

	//// Allocate memory for interleaver and de-interleaver tables ////
	int* Interleaver_table = (int*)malloc(Page_Size * Page_Size * sizeof(int)); // Interleaver table
	int* De_Interleaver_table = (int*)malloc(Page_Size * Page_Size * sizeof(int));  // De-interleaver table
	//// Allocate memory for interleaved data ////
	double** d_Interleaved = (double**)malloc(Page_Size * sizeof(double*)); // Interleaved double data
	int** n_Interleaved = (int**)malloc(Page_Size * sizeof(int*)); // Interleaved integer data
	for (i = 0; i < Page_Size; i++) {
		d_Interleaved[i] = (double*)calloc(Page_Size, sizeof(double));  // Allocate and initialize double interleaved data
		n_Interleaved[i] = (int*)calloc(Page_Size, sizeof(int));        // Allocate and initialize integer interleaved data
	}

	int** en_piece_sun_park;
	int** de_piece_sun_park;
	if (Encoder == 3)
	{
		en_piece_sun_park = (int**)malloc(64 * sizeof(int*));
		de_piece_sun_park = (int**)malloc(64 * sizeof(int*));
		for (i = 0; i < 64; i++)
		{
			en_piece_sun_park[i] = (int*)calloc(4, sizeof(int));
			de_piece_sun_park[i] = (int*)calloc(3, sizeof(int));
		}
		for (i = 0; i < 64; i++)
		{
			for (j = 0; j < 4; j++)
				en_piece_sun_park[i][j] = (i / int(pow(3.0, 3 - j))) % 3;
			for (j = 0; j < 3; j++)
				de_piece_sun_park[i][j] = (i >> ((2 - j) * 2)) % 4;
		}
	}

	FILE* OUTchi_24;


	int** piece_chi_24 = (int**)malloc(16 * sizeof(int*));
	for (i = 0; i < 16; i++)
		piece_chi_24[i] = (int*)calloc(4, sizeof(int));
	if (10 == Encoder)
	{
		OUTchi_24 = fopen("4ary24mc_16.txt", "r");
		for (i = 0; i < 16; i++)
			for (j = 0; j < 4; j++)
				fscanf(OUTchi_24, "%d", &piece_chi_24[i][j]);
		fclose(OUTchi_24);
	}
	//MakeInterleaverFile(Interleaver_table, De_Interleaver_table, Page_Size/3*Page_Size/2);

	snr = pow(10, (Mis_alignment_dB / 10.0));
	mis_sigma = sqrt(float(POWER / (snr)));

	HoloChannel(h, SigmaB, mis_sigma, SELL, Unit_Sell);
	PR_Target_set(Target_Coef, Target_Order);
	PRML_Convolution(trellis, Target_Coef[Target_Order / 2], Target_Order);

	if (Encoder == 6)
		PRML_Convolution_for_modul_gg_grad_23(trellis_modul);
	if (Encoder == 7)
	{
		PRML_Convolution_for_modul_gg_grad_23_prof(trellis_modul_prof);
		PRML_Convolution_for_modul_gg_grad_23_prof_method0(trellis_modul_prof_method0);
	}
	if (Encoder == 9)
	{
		PRML_Convolution_for_modul_gg_dmin2(trellis_modul_prof);
		PRML_Convolution_for_modul_gg_dmin2_method0(trellis_modul_prof_method0);
	}

	PRML_Convolution_gg_grad_23(trellis_A, Target_Coef[Target_Order / 2], Target_Order, 2, A_input);

	printf("Multilevel Hologram sigmaB %.2lf, Sell %d, unit_sell %d, Target_Order %d, Equalizer_Length %d\n", SigmaB, SELL, Unit_Sell, Target_Order, Equalizer_Length);
	printf("Quantization %d, page_size %d, block %d, total_bit %d\n", Quantization, Page_Size, BLOCK, long(Page_Size * Page_Size * BLOCK));
	printf("Adaptive_Length %d, Pridictor_Length %d, NoiseFilter %d, Adaptive %d, NF Adaptive %d\nDELTA %lf, ", Adaptive_Length, Pridictor_Length, NoiseFilter, Adaptive_yesornot, AdaptiveNF_yesornot, DELTA);
	printf("Detector is %s, ", s_detecter[DETECTOR]);
	if (Mis_const == 0)
		printf("Random Misalignment %lf db, \n", Mis_alignment_dB);
	else if (Mis_const == 1)
		printf("Misalignment (%1.2lf %1.2lf), \n", Mis_alignment_x, Mis_alignment_y);
	printf("Target is\n");
	for (i = 0; i < Target_Order; i++)
	{
		for (j = 0; j < Target_Order; j++)
			printf("%0.1lf ", Target_Coef[i][j]);
		printf("\n");
	}
	fprintf(fber, "Multilevel Hologram sigmaB %.2lf, Sell %d, unit_sell %d, Target_Order %d, Equalizer_Length %d\n", SigmaB, SELL, Unit_Sell, Target_Order, Equalizer_Length);
	fprintf(fber, "Quantization %d, page_size %d, block %d, total_bit %d\n", Quantization, Page_Size, BLOCK, long(Page_Size * Page_Size * BLOCK));
	fprintf(fber, "Adaptive_Length %d, Pridictor_Length %d, NoiseFilter %d, Adaptive %d, NF Adaptive %d\nDELTA %lf, ", Adaptive_Length, Pridictor_Length, NoiseFilter, Adaptive_yesornot, AdaptiveNF_yesornot, DELTA);
	fprintf(fber, "Detector is %s, ", s_detecter[DETECTOR]);
	if (Mis_const == 0)
		fprintf(fber, "Random Misalignment %lf db, \n", Mis_alignment_dB);
	else if (Mis_const == 1)
		fprintf(fber, "Misalignment (%1.2lf %1.2lf), \n", Mis_alignment_x, Mis_alignment_y);
	fprintf(fber, "Target is\n");
	for (i = 0; i < Target_Order; i++)
	{
		for (j = 0; j < Target_Order; j++)
			fprintf(fber, "%0.1lf ", Target_Coef[i][j]);
		fprintf(fber, "\n");
	}

	if (Encoder == 0)
	{
		printf("No modulation (random) code\n");
		fprintf(fber, "No modulation (random) code\n");
	}
	else if (Encoder == 1)
	{
		printf("Prof. Park 2/3 modulation code\n");
		fprintf(fber, "Prof. Park 2/3 modulation code\n");
	}
	else if (Encoder == 2)
	{
		printf("Kim's simple 3/4 modulation code\n");
		fprintf(fber, "Kim's simple 3/4 modulation code\n");
	}
	else if (Encoder == 3)
	{
		printf("Sunsuk & Gunhwan's 6/9 modulation code\n");
		fprintf(fber, "Sunsuk & Gunhwan's 6/9 modulation code\n");
	}
	else if (Encoder == 4)
	{
		printf("Kim odd even 2/3 modulation code\n");
		fprintf(fber, "Kim odd even 2/3 modulation code\n");
	}
	else if (Encoder == 5)
	{
		printf("Seungmin's 6/9 modulation code\n");
		fprintf(fber, "Seungmin's 6/9 modulation code\n");
	}
	else if (Encoder == 11)
	{
		printf("chi's 4/6 modulation code\n");
		//fprintf(fber, "Seungmin's 6/9 modulation code\n");
	}
	else if (Encoder == 6)
	{
		printf("gg_grad 2/3 modulation code\n");
		fprintf(fber, "gg_grad 2/3 modulation code\n");
	}
	else if (Encoder == 7)
	{
		printf("gg_grad 2/3 modulation code based on prof.park code\n");
		fprintf(fber, "gg_grad 2/3 modulation code based on prof.park code\n");
	}
	else if (Encoder == 8)
	{
		printf("gg_dmin2_code\n");
		fprintf(fber, "gg_dmin2_code\n");
	}
	else if (Encoder == 9)
	{
		printf("gg_dmin2_state_code based on modul method %d\n", modul_method);
		fprintf(fber, "gg_dmin2_state_code based on modul method %d\n", modul_method);
	}
	else if (Encoder == 10)
	{
		printf("Chi 24 proposal\n");
	}


	for (double SNR_DB = START_SNR - SNR_INTERVAL; SNR_DB < END_SNR;)
	{
		SNR_DB += SNR_INTERVAL;

		// Initialize variable
		count = 0, count_M_ary = 0, count_modul = 0, count_modul_M_ary = 0;
		block_count = 0, error = 0, error_eq = 0, error_eq_other = 0, error_M_ary = 0;
		error_viterbi[0] = 0, error_viterbi_M_ary[0] = 0;
		error_viterbi[1] = 0, error_viterbi_M_ary[1] = 0;
		error_demodul[0] = 0, error_demodul_M_ary[0] = 0;
		error_demodul[1] = 0, error_demodul_M_ary[1] = 0;

		for (i = 0; i < Equalizer_Length; i++)
			for (j = 0; j < Equalizer_Length; j++)
				Equalizer_Coef_x[i][j] = 0, Equalizer_Coef_y[i][j] = 0;

		snr = pow(10, (SNR_DB / 10.0)); // Công thức đã có, cứ vậy mà ốp
		sigma = float(POWER / (snr)); // Được dùng để mô phỏng nhiễu kênh

		while (block_count < Training_Block && DETECTOR != 3)
		{
			srand(0);
			if (Encoder == 0)
			{
				if (M_ary == 2)
				{
					for (i = 0; i < Data_Size; i++)
						for (j = 0; j < Data_Size; j++)
							inLDPC[i][j] = inLDPC_symbol[i][j] = rand() % 2, count++, count_M_ary++;
				}
				else if (M_ary == 4)
				{
					for (i = 0; i < Data_Size; i++)
						for (j = 0; j < Data_Size * 2; j++)
							inLDPC[i][j] = rand() % 2, count++;
					for (i = 0; i < Data_Size; i++)
						for (j = 0; j < Data_Size; j++)
							inLDPC_symbol[i][j] = (inLDPC[i][j * 2] << 1) + inLDPC[i][j * 2 + 1], count_M_ary++;
				}
			}



			else if (Encoder == 11)
			{
				if (M_ary != 4)
				{
					printf("\n\n mismatch encoder and M_ary\n You should select Encoder 3 and 4-ary.\n\n");
					exit(0);
				}
				for (i = 0; i < (Data_Size / 3) * (Data_Size / 2) * 4; i++)
					input_chi_46_4ary[i] = rand() % M_ary, count_modul += 2, count_modul_M_ary++;
				Encode_chiProposal46_4ary(inLDPC_symbol, input_chi_46_4ary, Page_Size);
				count_M_ary += Page_Size * Page_Size;
				count += Page_Size * Page_Size * 2;
			}

			for (i = 0; i < Page_Size + (Equalizer_Length / 2) * 2; i++)
				for (j = 0; j < Page_Size + (Equalizer_Length / 2) * 2; j++)
					if ((i > Equalizer_Length / 2 - 1 && i < Page_Size + Equalizer_Length / 2) && (j > Equalizer_Length / 2 - 1 && j < Page_Size + Equalizer_Length / 2))
						inCH[i][j] = inLDPC_symbol[i - Equalizer_Length / 2][j - Equalizer_Length / 2];
					else
						inCH[i][j] = 0;

			Hologram_ChannelOut(outCH, h, inCH, Page_Size + (Equalizer_Length / 2) * 2, SELL, sigma);
			//Hologram_EqualizerOut_x(outEQ_x, Equalizer_Coef_x, outCH, Page_Size, SELL, Equalizer_Length, 0, inCH, Target_Coef, Target_Order);
			Hologram_EqualizerOut_y(outEQ_y, Equalizer_Coef_y, outCH, Page_Size, SELL, Equalizer_Length, 0, inCH, Target_Coef, Target_Order);
			if (NoiseFilter == 1)
			{
				NPML_Pridictor_x(Pridictor_Coef_x, outEQ_x, Page_Size, SELL, Equalizer_Length, inCH, Target_Coef, Target_Order, Pridictor_Length);
				NPML_Pridictor_y(Pridictor_Coef_y, outEQ_y, Page_Size, SELL, Equalizer_Length, inCH, Target_Coef, Target_Order, Pridictor_Length);
			}
			block_count++;
		}
		block_count = 0, count = 0, count_M_ary = 0;
		count_modul = 0, count_modul_M_ary = 0;
		while (block_count < BLOCK)
		{
			if ((block_count + 1) % 50 == 0)
				printf(">");
			if ((block_count + 1) % 200 == 0)
				printf("%d", (block_count + 1) / 100);
			srand(block_count);
			if (Encoder == 0)
			{
				if (M_ary == 2)
				{
					for (i = 0; i < Data_Size; i++)
						for (j = 0; j < Data_Size; j++)
							inLDPC[i][j] = inLDPC_symbol[i][j] = rand() % 2, count++, count_M_ary++;
				}
				else if (M_ary == 4)
				{
					for (i = 0; i < Data_Size; i++)
						for (j = 0; j < Data_Size * 2; j++)
							inLDPC[i][j] = rand() % 2, count++;
					for (i = 0; i < Data_Size; i++)
						for (j = 0; j < Data_Size; j++)
							inLDPC_symbol[i][j] = (inLDPC[i][j * 2] << 1) + inLDPC[i][j * 2 + 1], count_M_ary++;
				}
			}


			else if (Encoder == 11)
			{
				if (M_ary != 4)
				{
					printf("\n\n mismatch encoder and M_ary\n You should select Encoder 3 and 4-ary.\n\n");
					exit(0);
				}
				for (i = 0; i < (Data_Size / 3) * (Data_Size / 2) * 4; i++)
					input_chi_46_4ary[i] = rand() % M_ary, 
					count_modul += 2, 
					count_modul_M_ary++;
				Encode_chiProposal46_4ary(inLDPC_symbol, input_chi_46_4ary, Page_Size);
				count_M_ary += Page_Size * Page_Size;
				count += Page_Size * Page_Size * 2;
			}

			if (file_see == 1)
			{
				sprintf(strH, "inLDPC_symbol_%1d.dat", block_count);
				OUT = fopen(strH, "w");
				for (i = 0; i < Page_Size; i++)
				{
					for (j = 0; j < Page_Size; j++)
						fprintf(OUT, "%d ", inLDPC_symbol[i][j]);
					fprintf(OUT, "\n");
				}
				fclose(OUT);
			}

			for (i = 0; i < Page_Size + (Equalizer_Length / 2) * 2; i++)
				for (j = 0; j < Page_Size + (Equalizer_Length / 2) * 2; j++)
					if ((i > Equalizer_Length / 2 - 1 && i < Page_Size + Equalizer_Length / 2) && (j > Equalizer_Length / 2 - 1 && j < Page_Size + Equalizer_Length / 2))
						inCH[i][j] = inLDPC_symbol[i - Equalizer_Length / 2][j - Equalizer_Length / 2];
					else
						inCH[i][j] = 0;

			if (file_see == 1)
			{
				OUT = fopen("inCH.dat", "w");
				for (i = 0; i < Page_Size + (Equalizer_Length / 2) * 2; i++)
				{
					for (j = 0; j < Page_Size + (Equalizer_Length / 2) * 2; j++)
						fprintf(OUT, "%d ", inCH[i][j]);
					fprintf(OUT, "\n");
				}
				fclose(OUT);
			}

			/////////////////////////////// channel out ///////////////////////////////////////////////
			HoloChannel(h, SigmaB, mis_sigma, SELL, Unit_Sell);
			Hologram_ChannelOut(outCH, h, inCH, Page_Size + (Equalizer_Length / 2) * 2, SELL, sigma);

			// Tính toán đầu ra sau khi nhận được kết quả nằm ở outCH
			for (i = 0; i < M_ary; i++)
				temp_count[i] = 0, temp_sum[i] = 0;
			for (i = Equalizer_Length / 2; i < Page_Size + (Equalizer_Length / 2); i++)
			{
				for (j = Equalizer_Length / 2; j < Page_Size + (Equalizer_Length / 2); j++)
					for (k = 0; k < M_ary; k++)
						if (inLDPC_symbol[i - Equalizer_Length / 2][j - Equalizer_Length / 2] == k)
							temp_sum[k] += outCH[i][j], temp_count[k]++;
			}
			for (k = 0; k < M_ary - 1; k++)
				chout_average[k] = ((temp_sum[k] / temp_count[k]) + (temp_sum[k + 1] / temp_count[k + 1])) / 2;
			//chout_average=(temp_sum_0+temp_sum_1)/(temp_count_0+temp_count_1);

			if (checkHist == 1)
			{
				for (k = 0; k < M_ary; k++)
				{
					sprintf(strH, "chout_%d.dat", k);
					CHOUT[k] = fopen(strH, "a");
				}
				for (i = 0; i < Page_Size; i++)
				{
					for (j = 0; j < Page_Size; j++)
					{
						for (k = 0; k < M_ary; k++)
						{
							if (inLDPC_symbol[i][j] == k)
							{
								fprintf(CHOUT[k], "%.10lf \n", outCH[i + Equalizer_Length / 2][j + Equalizer_Length / 2] / 10);
								break;
							}

						}
					}
				}
				for (k = 0; k < M_ary; k++)
					fclose(CHOUT[k]);
			}

			/////////////////////////////// equalization ///////////////////////////////////////////////
			//Hologram_EqualizerOut_x(outEQ_x, Equalizer_Coef_x, outCH, Page_Size, SELL, Equalizer_Length, 1, inCH, Target_Coef, Target_Order);
			Hologram_EqualizerOut_y(outEQ_y, Equalizer_Coef_y, outCH, Page_Size, SELL, Equalizer_Length, 1, inCH, Target_Coef, Target_Order);

			if (file_see == 1)
			{
				for (k = 0; k < M_ary; k++)
				{
					sprintf(strH, "eqout_%d.dat", k);
					CHOUT[k] = fopen(strH, "a");
				}
				for (i = 0; i < Page_Size; i++)
				{
					for (j = 0; j < Page_Size; j++)
					{
						for (k = 0; k < M_ary; k++)
						{
							if (inLDPC_symbol[i][j] == k)
							{
								fprintf(CHOUT[k], "%.10lf \n", outEQ_y[i][j]);
								break;
							}

						}
					}
				}
				for (k = 0; k < M_ary; k++)
					fclose(CHOUT[k]);
			}

			/////////////////////////////// detection ///////////////////////////////////////////////
			if (DETECTOR == 0)
			{
				//Hologram_SOVA_x(d_outViterbi_x, outEQ_x, EI, trellis, Pridictor_Coef_x, Pridictor_Length, Target_Order, Page_Size, depth);
				Hologram_SOVA_y(d_outViterbi_y, outEQ_y, EI, trellis, Pridictor_Coef_y, Pridictor_Length, Target_Order, Page_Size, depth);

			}
			else if (DETECTOR == 7)
			{
				Hologram_SOVA_x_for_gg_grad(d_outViterbi_x, outEQ_x, EI, trellis, trellis_A, Pridictor_Coef_x, Pridictor_Length, Target_Order, Page_Size, depth);
				ggg = 0;
			}

			if (file_see == 1)
			{
				OUT = fopen("d_outViterbi_x.dat", "w");
				for (i = 0; i < Page_Size; i++)
				{
					for (j = 0; j < Page_Size; j++)
						fprintf(OUT, "%lf ", d_outViterbi_x[i][j]);
					fprintf(OUT, "\n");
				}
				fclose(OUT);

				OUT = fopen("d_outViterbi_y.dat", "w");
				for (i = 0; i < Page_Size; i++)
				{
					for (j = 0; j < Page_Size; j++)
						fprintf(OUT, "%lf ", d_outViterbi_y[i][j]);
					fprintf(OUT, "\n");
				}
				fclose(OUT);
			}

			/////////////////////////////// demodulation //////////////////////////////////////////////	

			if (Encoder == 11)
			{
				/*Decode_chiProposal46_4ary(output_chi_46_4ary, d_outViterbi_x, Page_Size);
				for (i = 0; i < (Data_Size / 3)*(Data_Size / 2) * 4; i++)
				{
					if (input_chi_46_4ary[i] != output_chi_46_4ary[i])
						error_demodul_M_ary[0]++;
					if (input_chi_46_4ary[i] >> 1 != output_chi_46_4ary[i] >> 1)
						error_demodul[0]++;
					if (input_chi_46_4ary[i] % 2 != output_chi_46_4ary[i] % 2)
						error_demodul[0]++;
				}	*/
				Decode_chiProposal46_4ary(output_chi_46_4ary, d_outViterbi_y, Page_Size);
				for (i = 0; i < (Data_Size / 3) * (Data_Size / 2) * 4; i++)
				{
					if (input_chi_46_4ary[i] != output_chi_46_4ary[i])
						error_demodul_M_ary[1]++;
					if (input_chi_46_4ary[i] >> 1 != output_chi_46_4ary[i] >> 1)
						error_demodul[1]++;
					if (input_chi_46_4ary[i] % 2 != output_chi_46_4ary[i] % 2)
						error_demodul[1]++;
				}
			}


			/////////////////////////////// error check ///////////////////////////////////////////////
			if (M_ary == 4)
			{
				for (i = 0; i < Page_Size; i++)
				{
					for (j = 0; j < (Page_Size); j++)
					{
						ggg = 0;
						for (k = 0; k < M_ary - 1; k++)
						{
							if (outCH[i + Equalizer_Length / 2][j + Equalizer_Length / 2] < chout_average[k])
							{
								temp_count[0] = k;
								break;
							}
							else
								temp_count[0] = M_ary - 1;
						}
						if (inLDPC_symbol[i][j] != temp_count[0])
							error_M_ary++;
						if (inLDPC_symbol[i][j] >> 1 != temp_count[0] >> 1)
							error++;
						if (inLDPC_symbol[i][j] % 2 != temp_count[0] % 2)
							error++;

						//if (inLDPC_symbol[i][j] != (int)d_outViterbi_x[i][j])
						//	error_viterbi_M_ary[0]++;
						//if (inLDPC_symbol[i][j] >> 1 != (int)d_outViterbi_x[i][j] >> 1)
						//	error_viterbi[0]++;
						//if (inLDPC_symbol[i][j] % 2 != (int)d_outViterbi_x[i][j] % 2)
						//	error_viterbi[0]++;

						if (inLDPC_symbol[i][j] != (int)d_outViterbi_y[i][j])
							error_viterbi_M_ary[1]++;
						if (inLDPC_symbol[i][j] >> 1 != (int)d_outViterbi_y[i][j] >> 1)
							error_viterbi[1]++;
						if (inLDPC_symbol[i][j] % 2 != (int)d_outViterbi_y[i][j] % 2)
							error_viterbi[1]++;

					}

				}
			}

			block_count++;
		}
		printf("\n");
		//printf("count : %d\n count_M_ary : %d\n count_modul : %d\n count_modul_M_ary : %d\n\n",count, count_M_ary, count_modul, count_modul_M_ary);
		printf("snr_db : %.1lf\n", SNR_DB);
		//printf("snr_db : %.1lf\t\tMis snr_db : %.1lf\n", SNR_DB, Mis_alignment_dB);
		printf("channel error : %d\tBER : %.12lf\n", error, (double)error / count);
		printf("channel M_ary error : %d\tBER : %.12lf\n", error_M_ary, (double)error_M_ary / count_M_ary);
		printf("viterbi error x : %d\tBER : %.12lf\n", error_viterbi[0], (double)error_viterbi[0] / count);
		printf("viterbi M_ary error x : %d\tBER : %.12lf\n", error_viterbi_M_ary[0], (double)error_viterbi_M_ary[0] / count_M_ary);
		printf("viterbi error y : %d\tBER : %.12lf\n", error_viterbi[1], (double)error_viterbi[1] / count);
		printf("viterbi M_ary error y : %d\tBER : %.12lf\n", error_viterbi_M_ary[1], (double)error_viterbi_M_ary[1] / count_M_ary);
		if (Encoder != 0)
		{
			//printf("modulation error x : %d\tBER : %.12lf\n",error_demodul[0],(double)error_demodul[0]/count_modul);
			//printf("modulation M_ary error x : %d\tBER : %.12lf\n",error_demodul_M_ary[0],(double)error_demodul_M_ary[0]/count_modul_M_ary);
			printf("modulation error y : %d\tBER : %.12lf\n", error_demodul[1], (double)error_demodul[1] / count_modul);
			printf("modulation M_ary error y : %d\tBER : %.12lf\n", error_demodul_M_ary[1], (double)error_demodul_M_ary[1] / count_modul_M_ary);
		}

		//fprintf(fber, "\n");
		//fprintf(fber, "snr_db : %.1lf\t\tMis snr_db : %.1lf\n", SNR_DB, Mis_alignment_dB);
		fprintf(fber_outCH, "%d\t%.12lf\n", int(SNR_DB / SNR_INTERVAL), (double)error / count);
		fprintf(fber_outCH_Mary, "%d\t%.12lf\n", int(SNR_DB / SNR_INTERVAL), (double)error_M_ary / count_M_ary);
		//fprintf(fber,"viterbi error x : %d\tBER : %.12lf\n",error_viterbi[0],(double)error_viterbi[0]/count);
		//fprintf(fber,"viterbi M_ary error x : %d\tBER : %.12lf\n",error_viterbi_M_ary[0],(double)error_viterbi_M_ary[0]/count_M_ary);
		fprintf(fber_outSOVA, "%d\t%.12lf\n", int(SNR_DB / SNR_INTERVAL), (double)error_viterbi[1] / count);
		fprintf(fber_outSOVA_Mary, "%d\t%.12lf\n", int(SNR_DB / SNR_INTERVAL), (double)error_viterbi_M_ary[1] / count_M_ary);
		if (Encoder != 0)
		{
			//fprintf(fber,"modulation error x : %d\tBER : %.12lf\n",error_demodul[0],(double)error_demodul[0]/count);
			//fprintf(fber,"modulation M_ary error x : %d\tBER : %.12lf\n",error_demodul_M_ary[0],(double)error_demodul_M_ary[0]/count_modul_M_ary);
			fprintf(fber_outDEM, "%d\t%.12lf\n", int(SNR_DB / SNR_INTERVAL), (double)error_demodul[1] / count_modul);
			fprintf(fber_outDEM_Mary, "%d\t%.12lf\n", int(SNR_DB / SNR_INTERVAL), (double)error_demodul_M_ary[1] / count_modul_M_ary);
		}
		if (DRAWING_BER_CURVE == 0)
		{
			fprintf(BRIEF_RESULT, "%s_channel(%d)=%.12lf;\n", s_encoder[Encoder], int(SNR_DB / SNR_INTERVAL), (double)error / count);
			fprintf(BRIEF_RESULT, "%s_channel_M_ary(%d)=%.12lf;\n", s_encoder[Encoder], int(SNR_DB / SNR_INTERVAL), (double)error_M_ary / count_M_ary);
			fprintf(BRIEF_RESULT, "%%data(x)\n");
			//fprintf(BRIEF_RESULT,"%s_%s(%d,:)=%.12lf;\n",s_encoder[Encoder],s_detecter[DETECTOR],int(SNR_DB/SNR_INTERVAL),(double)error_viterbi[0]/count);
			//fprintf(BRIEF_RESULT,"%s_%s_M_ary(%d,:)=%.12lf;\n",s_encoder[Encoder],s_detecter[DETECTOR],int(SNR_DB/SNR_INTERVAL),(double)error_viterbi_M_ary[0]/count_M_ary);
			fprintf(BRIEF_RESULT, "%s_%s(%d,:)=%.12lf;\n", s_encoder[Encoder], s_detecter[DETECTOR], int(SNR_DB / SNR_INTERVAL), (double)error_viterbi[1] / count);
			fprintf(BRIEF_RESULT, "%s_%s_M_ary(%d,:)=%.12lf;\n", s_encoder[Encoder], s_detecter[DETECTOR], int(SNR_DB / SNR_INTERVAL), (double)error_viterbi_M_ary[1] / count_M_ary);
			if (Encoder != 0)
			{
				//fprintf(BRIEF_RESULT,"%s_%s_demodul(%d,:)=%.12lf;\n",s_encoder[Encoder],s_detecter[DETECTOR],int(SNR_DB/SNR_INTERVAL),(double)error_demodul[0]/count_modul);
				//fprintf(BRIEF_RESULT,"%s_%s_demodul_M_ary(%d,:)=%.12lf;\n",s_encoder[Encoder],s_detecter[DETECTOR],int(SNR_DB/SNR_INTERVAL),(double)error_demodul_M_ary[0]/count_modul_M_ary);
				fprintf(BRIEF_RESULT, "%s_%s_demodul(%d,:)=%.12lf;\n", s_encoder[Encoder], s_detecter[DETECTOR], int(SNR_DB / SNR_INTERVAL), (double)error_demodul[1] / count_modul);
				fprintf(BRIEF_RESULT, "%s_%s_demodul_M_ary(%d,:)=%.12lf;\n", s_encoder[Encoder], s_detecter[DETECTOR], int(SNR_DB / SNR_INTERVAL), (double)error_demodul_M_ary[1] / count_modul_M_ary);
			}
		}
		/*
		if(DRAWING_BER_CURVE==1)
		{
			fprintf(BRIEF_RESULT,"%s_channel(%d)=%.12lf;\n",s_encoder[Encoder],int(SNR_DB/SNR_INTERVAL),(double)error/count);
			fprintf(BRIEF_RESULT,"%s_channel_M_ary(%d)=%.12lf;\n",s_encoder[Encoder],int(SNR_DB/SNR_INTERVAL),(double)error_M_ary/count_M_ary);
			fprintf(BRIEF_RESULT,"%%data(x,y)\n");
			fprintf(BRIEF_RESULT,"%s_%s(%d,:)=[%.12lf %.12lf];\n",s_encoder[Encoder],s_detecter[DETECTOR],int(SNR_DB/SNR_INTERVAL),(double)error_viterbi[0]/count,(double)error_viterbi[1]/count);
			fprintf(BRIEF_RESULT,"%s_%s_M_ary(%d,:)=[%.12lf %.12lf];\n",s_encoder[Encoder],s_detecter[DETECTOR],int(SNR_DB/SNR_INTERVAL),(double)error_viterbi_M_ary[0]/count_M_ary,(double)error_viterbi_M_ary[1]/count_M_ary);
			if(Encoder != 0)
			{
				fprintf(BRIEF_RESULT,"%s_%s_demodul(%d,:)=[%.12lf %.12lf];\n",s_encoder[Encoder],s_detecter[DETECTOR],int(SNR_DB/SNR_INTERVAL),(double)error_demodul[0]/count_modul,(double)error_demodul[1]/count_modul);
				fprintf(BRIEF_RESULT,"%s_%s_demodul_M_ary(%d,:)=[%.12lf %.12lf];\n",s_encoder[Encoder],s_detecter[DETECTOR],int(SNR_DB/SNR_INTERVAL),(double)error_demodul_M_ary[0]/count_modul_M_ary,(double)error_demodul_M_ary[1]/count_modul_M_ary);
			}
		}
		*/
	}


	c_time = clock() - c_time; // c_time Check End 
	//fprintf(fber, "%d½Ã %dºÐ %dÃÊÀÔ´Ï´Ù. \n", int(c_time / CLOCKS_PER_SEC) / 3600, (int(c_time / CLOCKS_PER_SEC) % 3600) / 60, int(c_time / CLOCKS_PER_SEC) % 60);
	printf("%dhour %dmin %dsec. \n", int(c_time / CLOCKS_PER_SEC) / 3600, (int(c_time / CLOCKS_PER_SEC) % 3600) / 60, int(c_time / CLOCKS_PER_SEC) % 60);

	getchar();

	fcloseall();

	for (i = 0; i < SELL; i++)
		free(h[i]);
	free(h);

	for (i = 0; i < Target_Order; i++)
		free(Target_Coef[i]);
	free(Target_Coef);

	for (i = 0; i < Page_Size; i++)
		free(inLDPC[i]), free(EI[i]), free(outLDPC[i]);
	free(inLDPC), free(EI), free(outLDPC);

	for (i = 0; i < pow((double)M_ary, Target_Order - 1); i++)
		free(trellis[i]);
	free(trellis);

	for (i = 0; i < Equalizer_Length; i++)
		free(Equalizer_Coef_x[i]), free(Equalizer_Coef_y[i]);
	free(Equalizer_Coef_x), free(Equalizer_Coef_y);

	for (i = 0; i < Page_Size + (Equalizer_Length / 2) * 2; i++)
		free(inCH[i]), free(outCH[i]);
	free(inCH), free(outCH);

	for (i = 0; i < Page_Size; i++)
		free(outEQ_x[i]), free(outEQ_y[i]), free(n_outViterbi[i]), free(d_outViterbi[i]), free(d_outViterbi_x[i]), free(d_outViterbi_y[i]);
	free(d_outViterbi), free(n_outViterbi), free(outEQ_x), free(outEQ_y);
	free(d_outViterbi_x), free(d_outViterbi_y);


	free(input_park_2_3_4ary), free(output_park_2_3_4ary);

	for (i = 0; i < M_ary - 2; i++)
		free(trellis_A[i]);
	free(trellis_A);

	for (i = 0; i < (int)pow((double)M_ary, 2); i++)
		for (j = 0; j < (int)pow((double)M_ary, 2); j++)
			free(trellis_modul[i][j]);
	for (i = 0; i < (int)pow((double)M_ary, 2); i++)
		free(trellis_modul[i]);
	free(trellis_modul);

	for (i = 0; i < (int)pow((double)M_ary, 2); i++)
		for (j = 0; j < (int)pow((double)M_ary, 2); j++)
			free(trellis_modul_prof[i][j]);
	for (i = 0; i < (int)pow((double)M_ary, 2); i++)
		free(trellis_modul_prof[i]);
	free(trellis_modul_prof);

	for (i = 0; i < total_states_method0; i++)
		for (j = 0; j < total_states_method0; j++)
			for (k = 0; k < 2; k++)
				free(trellis_modul_prof_method0[i][j][k]);
	for (i = 0; i < total_states_method0; i++)
		for (j = 0; j < total_states_method0; j++)
			free(trellis_modul_prof_method0[i][j]);
	for (i = 0; i < total_states_method0; i++)
		free(trellis_modul_prof_method0[i]);
	free(trellis_modul_prof_method0);




}

// 1
void HoloChannel(double** h, double SigmaB, double Mis_alignment_sigma, int SELL, int Unit_Sell)
{
	int i, j, l, m;

	double mis_length = gaussian(0.0, Mis_alignment_sigma, &k_noise);
	double mis_radian = gaussian(0.0, Mis_alignment_sigma, &k_noise);
	double mis_x = mis_length * cos(mis_radian);
	double mis_y = mis_length * sin(mis_radian);

	if (Mis_const == 1)
		mis_x = Mis_alignment_x, mis_y = Mis_alignment_y;;


	double** temp_h = (double**)malloc(SELL * Unit_Sell * sizeof(double*));
	for (i = 0; i < SELL * Unit_Sell; i++)
		temp_h[i] = (double*)calloc(SELL * Unit_Sell, sizeof(double));

	for (i = -SELL / 2; i <= SELL / 2; i++)
	{
		for (j = -SELL / 2; j <= SELL / 2; j++)
		{
			for (l = -Unit_Sell / 2; l <= Unit_Sell / 2; l++)
			{
				for (m = -Unit_Sell / 2; m <= Unit_Sell / 2; m++)
				{
					if ((i + mis_x) * Unit_Sell + l == 0 && (j + mis_y) * Unit_Sell + m == 0)
						temp_h[i * Unit_Sell + l + (Unit_Sell * (SELL / 2) + Unit_Sell / 2)][j * Unit_Sell + m + (Unit_Sell * (SELL / 2) + Unit_Sell / 2)] = 1 / (SigmaB * SigmaB);
					else if ((i + mis_x) * Unit_Sell + l == 0 && (j + mis_y) * Unit_Sell + m != 0)
						temp_h[i * Unit_Sell + l + (Unit_Sell * (SELL / 2) + Unit_Sell / 2)][j * Unit_Sell + m + (Unit_Sell * (SELL / 2) + Unit_Sell / 2)] = 1 / (SigmaB * SigmaB) * pow((sin(PI * (j + m / float(Unit_Sell) + mis_y) / SigmaB) / (PI * (j + m / float(Unit_Sell) + mis_y) / SigmaB)), 2);
					else if ((i + mis_x) * Unit_Sell + l != 0 && (j + mis_y) * Unit_Sell + m == 0)
						temp_h[i * Unit_Sell + l + (Unit_Sell * (SELL / 2) + Unit_Sell / 2)][j * Unit_Sell + m + (Unit_Sell * (SELL / 2) + Unit_Sell / 2)] = 1 / (SigmaB * SigmaB) * pow((sin(PI * (i + l / float(Unit_Sell) + mis_x) / SigmaB) / (PI * (i + l / float(Unit_Sell) + mis_x) / SigmaB)), 2);
					else
						temp_h[i * Unit_Sell + l + (Unit_Sell * (SELL / 2) + Unit_Sell / 2)][j * Unit_Sell + m + (Unit_Sell * (SELL / 2) + Unit_Sell / 2)] = 1 / (SigmaB * SigmaB) * pow((sin(PI * (i + l / float(Unit_Sell) + mis_x) / SigmaB) / (PI * (i + l / float(Unit_Sell) + mis_x) / SigmaB)) * (sin(PI * (j + m / float(Unit_Sell) + mis_y) / SigmaB) / (PI * (j + m / float(Unit_Sell) + mis_y) / SigmaB)), 2);
				}
			}
		}
	}
	for (i = 0; i < SELL; i++)
		for (j = 0; j < SELL; j++)
		{
			h[i][j] = 0;
			for (l = 0; l < Unit_Sell; l++)
				for (m = 0; m < Unit_Sell; m++)
					h[i][j] += temp_h[i * Unit_Sell + l][j * Unit_Sell + m];
		}

	for (i = 0; i < SELL * Unit_Sell; i++)
		free(temp_h[i]);
	free(temp_h);
}
// 3
void PR_Target_set(double** Target_Coef, int target_order)											//// Setting the coefficients of the PR Target according to the channel model.
{
	// Ma trận 3x3
	// 010
	// 131
	// 010
	if (target_order == 3)
	{
		Target_Coef[0][0] = 0, Target_Coef[0][1] = 1, Target_Coef[0][2] = 0;
		Target_Coef[1][0] = 1, Target_Coef[1][1] = 3, Target_Coef[1][2] = 1;
		//		if(DETECTOR==1)
		//			Target_Coef[1][0]=1, Target_Coef[1][1]=3, Target_Coef[1][2]=1;
		//		else if(DETECTOR==5)
		//			Target_Coef[1][0]=1, Target_Coef[1][1]=4, Target_Coef[1][2]=1;
		//		else
		//			Target_Coef[1][0]=1, Target_Coef[1][1]=4, Target_Coef[1][2]=1;
		//		Target_Coef[1][0]=1, Target_Coef[1][1]=int(Mis_alignment_dB), Target_Coef[1][2]=1;
		Target_Coef[2][0] = 0, Target_Coef[2][1] = 1, Target_Coef[2][2] = 0;
	}
	else if (target_order == 5)
	{
		Target_Coef[0][0] = 0, Target_Coef[0][1] = 0, Target_Coef[0][2] = 1, Target_Coef[0][3] = 0, Target_Coef[0][4] = 0;
		Target_Coef[1][0] = 0, Target_Coef[1][1] = 1, Target_Coef[1][2] = 2, Target_Coef[1][3] = 1, Target_Coef[1][4] = 0;
		Target_Coef[2][0] = 1, Target_Coef[2][1] = 2, Target_Coef[2][2] = 8, Target_Coef[2][3] = 2, Target_Coef[2][4] = 1;
		Target_Coef[3][0] = 0, Target_Coef[3][1] = 1, Target_Coef[3][2] = 2, Target_Coef[3][3] = 1, Target_Coef[3][4] = 0;
		Target_Coef[4][0] = 0, Target_Coef[4][1] = 0, Target_Coef[4][2] = 1, Target_Coef[4][3] = 0, Target_Coef[4][4] = 0;
	}
	else if (target_order == 7)
	{
		Target_Coef[0][0] = 0, Target_Coef[0][1] = 0, Target_Coef[0][2] = 0, Target_Coef[0][3] = 0, Target_Coef[0][4] = 0, Target_Coef[0][5] = 0, Target_Coef[0][6] = 0;
		Target_Coef[1][0] = 0, Target_Coef[1][1] = 0, Target_Coef[1][2] = 0, Target_Coef[1][3] = 0, Target_Coef[1][4] = 0, Target_Coef[1][5] = 0, Target_Coef[1][6] = 0;
		Target_Coef[2][0] = 0, Target_Coef[2][1] = 0, Target_Coef[2][2] = 0, Target_Coef[2][3] = 0, Target_Coef[2][4] = 0, Target_Coef[2][5] = 0, Target_Coef[2][6] = 0;
		Target_Coef[3][0] = 0, Target_Coef[3][1] = 0, Target_Coef[3][2] = 0, Target_Coef[3][3] = 0, Target_Coef[3][4] = 0, Target_Coef[3][5] = 0, Target_Coef[3][6] = 0;
		Target_Coef[4][0] = 0, Target_Coef[4][1] = 0, Target_Coef[4][2] = 0, Target_Coef[4][3] = 0, Target_Coef[4][4] = 0, Target_Coef[4][5] = 0, Target_Coef[4][6] = 0;
		Target_Coef[5][0] = 0, Target_Coef[5][1] = 0, Target_Coef[5][2] = 0, Target_Coef[5][3] = 0, Target_Coef[5][4] = 0, Target_Coef[5][5] = 0, Target_Coef[5][6] = 0;
		Target_Coef[6][0] = 0, Target_Coef[6][1] = 0, Target_Coef[6][2] = 0, Target_Coef[6][3] = 0, Target_Coef[6][4] = 0, Target_Coef[6][5] = 0, Target_Coef[6][6] = 0;
	}
}
void PRML_Convolution(double** trellis, double* Target_Coef, int Target_Order)
{
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//// trellis[i][j]	i=state, j= 0,1 => reference_bit(0 up, 1 down) 2,3 => reference_value(2 up, 3 down) ////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Target_Order--;
	int i, j, k;					//// Temporary variables used during the for loop
	int total_states = (int)pow((double)M_ary, Target_Order);


	int* temp_trel_a = (int*)calloc(Target_Order, sizeof(int));

	//// setting reference_bit ////
	for (i = 0; i < total_states; i++)
	{
		for (j = 0; j < M_ary; j++)
			trellis[i][j] = i % M_ary;
	}
	///////////////////////////////

	for (i = 0; i < Target_Order; i++)
		temp_trel_a[i] = 0;										// Initialize to 0
	for (i = 0; i < total_states; i++)
	{
		for (k = 0; k < Target_Order; k++)
			temp_trel_a[k] = (i / (int)pow((double)M_ary, k)) % M_ary;
		for (k = 0; k < M_ary; k++)
		{
			trellis[i][M_ary + k] = 0;
			for (j = 0; j < Target_Order; j++)
				trellis[i][M_ary + k] += Target_Coef[j] * temp_trel_a[j];
			trellis[i][M_ary + k] += (Target_Coef[Target_Order] * k);
			trellis[i][M_ary + k] /= M_ary - 1;
		}
	}


	free(temp_trel_a);
}
// 2
void Hologram_ChannelOut(double** outCH, double** h, int** inCH, int Page_Size, int SELL, double sigma)
{
	int i, j, x, y;
	double out;

	double sigma_AWGN = sqrt(sigma);

	for (x = 0; x < Page_Size; x++)
	{
		for (y = 0; y < Page_Size; y++)
		{
			out = 0;
			for (i = -SELL / 2; i <= SELL / 2; i++)
			{
				for (j = -SELL / 2; j <= SELL / 2; j++)
				{
					if (x + i >= 0 && y + j >= 0 && x + i < Page_Size && y + j < Page_Size)
					{
						out += inCH[x + i][y + j] * h[i + SELL / 2][j + SELL / 2];
					}
				}
			}
			if (noise_free == 0)
				outCH[x][y] = out / Unit_Sell / (M_ary - 1);
			else if (noise_free == 1)
				outCH[x][y] = out / Unit_Sell / (M_ary - 1) + gaussian(0.0, sigma_AWGN, &k_noise);
		}
	}
}
// 4
void Hologram_EqualizerOut_x(double** outEQ, double** eq_coef, double** inEQ, int Page_Size, int SELL, int Equalizer_Length, int Training, int** inCH, double** Target_Coef, int Target_Order)
{
	int i, j, l, m;

	double Target_Value;

	//FILE *ERROR=fopen("eqout_error.dat","a");
	for (i = 0; i < Page_Size; i++)
	{
		for (j = 0; j < Page_Size; j++)
		{
			outEQ[i][j] = 0;
			for (l = 0; l < Equalizer_Length; l++)
				//l=Equalizer_Length/2;
				for (m = 0; m < Equalizer_Length; m++)
					outEQ[i][j] += inEQ[i + l][j + m] * eq_coef[l][m];
			Target_Value = 0;
			if (Training == 0 && i > Equalizer_Length / 2 && j > Equalizer_Length / 2)
			{
				l = 0;
				//for(l=-Target_Order/2 ; l<=Target_Order/2 ; l++)
				for (m = -Target_Order / 2; m <= Target_Order / 2; m++)
					//						Target_Value += (inCH[i+Equalizer_Length/2+l][j+Equalizer_Length/2+m]==0 ? -1 : 1)*Target_Coef[l+Target_Order/2][m+Target_Order/2];
					Target_Value += (double)inCH[i + Equalizer_Length / 2 + l][j + Equalizer_Length / 2 + m] / (M_ary - 1) * Target_Coef[l + Target_Order / 2][m + Target_Order / 2];
				//Target_Value += (double)inCH[i+Equalizer_Length/2+l][j+Equalizer_Length/2+m]*Target_Coef[l+Target_Order/2][m+Target_Order/2];

				for (l = 0; l < Equalizer_Length; l++)
					//l=Equalizer_Length/2;
					for (m = 0; m < Equalizer_Length; m++)
						eq_coef[l][m] += inEQ[i + l][j + m] * (Target_Value - outEQ[i][j]) * DELTA;
			}
			//if(Training == 0)
			//	fprintf(ERROR,"%lf\n",Target_Value-outEQ[i][j]);			
		}
	}
	//fclose(ERROR);
}
void Hologram_EqualizerOut_y(double** outEQ, double** eq_coef, double** inEQ, int Page_Size, int SELL, int Equalizer_Length, int Training, int** inCH, double** Target_Coef, int Target_Order)
{
	int i, j, l, m;

	double Target_Value;

	//	FILE *ERROR=fopen("eqout_error.dat","w");
	for (j = 0; j < Page_Size; j++)
	{
		for (i = 0; i < Page_Size; i++)
		{
			outEQ[i][j] = 0;
			for (l = 0; l < Equalizer_Length; l++)
				//l=Equalizer_Length/2;
				for (m = 0; m < Equalizer_Length; m++)
					outEQ[i][j] += inEQ[i + l][j + m] * eq_coef[l][m];
			Target_Value = 0;
			if (Training == 0 && i > Equalizer_Length / 2 && j > Equalizer_Length / 2)
			{
				m = 0;
				for (l = -Target_Order / 2; l <= Target_Order / 2; l++)
					//	for(m=-Target_Order/2 ; m<=Target_Order/2 ; m++)
	//						Target_Value += (inCH[i+Equalizer_Length/2+l][j+Equalizer_Length/2+m]==0 ? -1 : 1)*Target_Coef[l+Target_Order/2][m+Target_Order/2];
					Target_Value += (double)inCH[i + Equalizer_Length / 2 + l][j + Equalizer_Length / 2 + m] / (M_ary - 1) * Target_Coef[l + Target_Order / 2][m + Target_Order / 2];

				for (l = 0; l < Equalizer_Length; l++)
					//l=Equalizer_Length/2;
					for (m = 0; m < Equalizer_Length; m++)
						eq_coef[l][m] += inEQ[i + l][j + m] * (Target_Value - outEQ[i][j]) * DELTA;
			}
			//			if(Training == 1)
			//				fprintf(ERROR,"%lf\n",Target_Value-outEQ[i][j]);			
		}
	}
}
// 5
void NPML_Pridictor_x(double* Pridictor_Coef_x, double** inEQ, int Page_Size, int SELL, int Equalizer_Length, int** inCH, double** Target_Coef, int Target_Order, int Pridictor_Length)
{
	int i, j, l, m;

	double Target_Value;
	double prierror, pout;

	//	int count=0;
	//	FILE *OUT, *OUT1;

	//	OUT=fopen("npin.dat","a");
	//	OUT1=fopen("nperror.dat","a");

	double* Pridictor_Memory = (double*)calloc((Pridictor_Length + 1), sizeof(double));

	for (i = 0; i < Page_Size; i++)
	{
		for (j = 0; j < Page_Size; j++)
		{
			pout = 0.0;
			for (l = 0; l < Pridictor_Length + 1; l++)
				pout += Pridictor_Memory[l] * Pridictor_Coef_x[l];

			Target_Value = 0;
			for (l = -Target_Order / 2; l < Target_Order / 2; l++)
				for (m = -Target_Order / 2; m < Target_Order / 2; m++)
					Target_Value += inCH[i + Equalizer_Length / 2 + l][j + Equalizer_Length / 2 + m] * Target_Coef[l + Target_Order / 2][m + Target_Order / 2];

			prierror = inEQ[i][j] - Target_Value - pout;

			//			if(count>100000 && count<500000)
			//			{
			//				fprintf(OUT,"%lf\n",inEQ[i][j] - Target_Value);
			//				fprintf(OUT1,"%lf\n",prierror);
			//			}
			//			count++;


			for (l = 0; l < Pridictor_Length + 1; l++)
				Pridictor_Coef_x[l] += DELTA * prierror * Pridictor_Memory[l];

			///////////////////	 Pridictor_Memory Buffer   ////////////////////////
			for (l = Pridictor_Length; l > 0; l--)
				Pridictor_Memory[l] = Pridictor_Memory[l - 1];
			Pridictor_Memory[0] = inEQ[i][j] - Target_Value;
		}
	}
	free(Pridictor_Memory);
}
void NPML_Pridictor_y(double* Pridictor_Coef_y, double** inEQ, int Page_Size, int SELL, int Equalizer_Length, int** inCH, double** Target_Coef, int Target_Order, int Pridictor_Length)
{
	int i, j, l, m;

	double Target_Value;
	double prierror, pout;

	double* Pridictor_Memory = (double*)calloc((Pridictor_Length + 1), sizeof(double));

	for (j = 0; j < Page_Size; j++)
	{
		for (i = 0; i < Page_Size; i++)
		{
			pout = 0.0;
			for (l = 0; l < Pridictor_Length + 1; l++)
				pout += Pridictor_Memory[l] * Pridictor_Coef_y[l];

			Target_Value = 0;
			for (l = -Target_Order / 2; l < Target_Order / 2; l++)
				for (m = -Target_Order / 2; m < Target_Order / 2; m++)
					Target_Value += inCH[i + Equalizer_Length / 2 + l][j + Equalizer_Length / 2 + m] * Target_Coef[l + Target_Order / 2][m + Target_Order / 2];

			prierror = inEQ[i][j] - Target_Value - pout;

			for (l = 0; l < Pridictor_Length + 1; l++)
				Pridictor_Coef_y[l] += DELTA * prierror * Pridictor_Memory[l];

			///////////////////	 Pridictor_Memory Buffer   ////////////////////////
			for (l = Pridictor_Length; l > 0; l--)
				Pridictor_Memory[l] = Pridictor_Memory[l - 1];
			Pridictor_Memory[0] = inEQ[i][j] - Target_Value;
		}
	}
	free(Pridictor_Memory);
}
// 6  2D soft-output Viterbi algorithm (SOVA) as the 2D detection for the HDS systems 
void Hologram_SOVA_x(double** viterbi_out, double** eq_out, double** EI, double** trellis, double* Pridictor_Coef, int Pridictor_Length, int Target_Order, int Page_Size, int depth)
{
	Target_Order--;
	int i, j, k, x, y;				//// Temporary variables used during the for loop
	int Viterbi_temp;
	int vout = 0;
	double pri_error[M_ary];
	int over;
	double PM[M_ary];
	double alpha = 1.0;
	int error_viterbi = 0;
	int scale = 30;
	int iteration = 3;
	int total_states = (int)pow((double)M_ary, Target_Order);

	int ggg = 0;


	double eqout_total = 0, eqout_average;
	for (i = 0; i < Page_Size; i++)
		for (j = 0; j < Page_Size; j++)
			eqout_total += eq_out[i][j];
	eqout_average = eqout_total / (Page_Size * Page_Size);

	double** differ = (double**)malloc(total_states * sizeof(double*));
	for (i = 0; i < total_states; i++)
		differ[i] = (double*)calloc(M_ary, sizeof(double));

	for (i = 0; i < Page_Size; i++)
		for (j = 0; j < Page_Size; j++)
			viterbi_out[i][j] = 10000;

	unsigned long int*** Viterbi_State = (unsigned long int***)malloc(total_states * sizeof(unsigned long int**));
	for (i = 0; i < total_states; i++)
		Viterbi_State[i] = (unsigned long int**)malloc(depth * sizeof(unsigned long int*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			Viterbi_State[i][j] = (unsigned long int*)calloc(2, sizeof(unsigned long int));

	double** Viterbi_Path = (double**)malloc(total_states * sizeof(double*));
	for (i = 0; i < total_states; i++)
		Viterbi_Path[i] = (double*)calloc(2, sizeof(double));

	double*** Viterbi_whereState = (double***)malloc(total_states * sizeof(double**));
	for (i = 0; i < total_states; i++)
		Viterbi_whereState[i] = (double**)malloc(depth * sizeof(double*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			Viterbi_whereState[i][j] = (double*)calloc(2, sizeof(double));

	double*** Viterbi_error_Memory = (double***)malloc(total_states * sizeof(double**));
	for (i = 0; i < total_states; i++)
		Viterbi_error_Memory[i] = (double**)malloc((Pridictor_Length + 2) * sizeof(double*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < Pridictor_Length + 2; j++)
			Viterbi_error_Memory[i][j] = (double*)calloc(2, sizeof(double));

	for (x = 0; x < Page_Size; x++)
	{
		for (i = 1; i < total_states; i++)
			Viterbi_Path[i][1] = 10000, Viterbi_State[i][0][0] = 0, Viterbi_State[i][0][1] = 0, Viterbi_State[i][1][0] = 0, Viterbi_State[i][1][1] = 0;
		Viterbi_Path[0][1] = 0, Viterbi_State[0][0][0] = 0, Viterbi_State[0][0][1] = 0, Viterbi_State[0][1][0] = 0, Viterbi_State[0][1][1] = 0;

		for (i = 0; i < total_states; i++)
			for (j = 0; j < Pridictor_Length + 2; j++)
				Viterbi_error_Memory[i][j][0] = 0.0, Viterbi_error_Memory[i][j][1] = 0.0;

		for (y = 0; y < Page_Size + 1; y++)
		{

			ggg = 0;
			for (i = 0; i < total_states; i++)
			{
				// calculate the predicted error
				pri_error[0] = 0.0, pri_error[1] = 0.0;
				for (j = 1; j < Pridictor_Length + 2; j++)
				{
					for (k = 0; k < M_ary; k++)
						pri_error[k] = Viterbi_error_Memory[i / M_ary + k * M_ary][j][1] * Pridictor_Coef[j - 1];
				}
				// calculate each path metric of a states
				if (NoiseFilter == 0)	for (k = 0; k < M_ary; k++) pri_error[k] = 0.0;
				if (y < Page_Size + 1 && y>0)
				{
					for (k = 0; k < M_ary; k++)
						PM[k] = Viterbi_Path[i / M_ary + k * M_ary][1] + pow(eq_out[x][y - 1] - trellis[i][M_ary + k] - pri_error[0], 2);// - (EI[x][y-1]/2)*(trellis[i][k]*2-1);
				}
				else
				{
					for (k = 0; k < M_ary; k++)
						PM[k] = Viterbi_Path[i / M_ary + k * M_ary][1];// +pow(eq_out[x][y - 1] - trellis[i][M_ary + k] - pri_error[0], 2);// - (EI[x][y-1]/2)*(trellis[i][k]*2-1);
				}

				// find minimum path of each state				
				Viterbi_temp = 0;
				if (gg_pr_viterbi == 1)
				{
					for (k = 0; k < M_ary; k++)
					{
						if (k == (M_ary - 1))
						{
							if ((i > 3) && (i < 13))
							{
								if (PM[Viterbi_temp] > PM[k])
									Viterbi_temp = k;
							}
						}
						else
						{
							if (PM[Viterbi_temp] > PM[k])
								Viterbi_temp = k;
						}
					}
				}
				else
				{
					for (k = 0; k < M_ary; k++)
						if (PM[Viterbi_temp] > PM[k])
							Viterbi_temp = k;
				}


				// save differents between survivor path and competition path
				for (k = 0; k < M_ary; k++)
					differ[i][k] = PM[k] - PM[Viterbi_temp];

				// save the survivor path of each state
				Viterbi_Path[i][0] = PM[Viterbi_temp];
				for (j = 1; j < depth; j++)
					Viterbi_State[i][j][0] = Viterbi_State[i / M_ary + Viterbi_temp * M_ary][j][1];
				Viterbi_State[i][0][0] = (unsigned long)trellis[i][Viterbi_temp];
				for (j = 1; j < depth; j++)
					Viterbi_whereState[i][j][0] = Viterbi_whereState[i / M_ary + Viterbi_temp * M_ary][j][1];
				Viterbi_whereState[i][0][0] = Viterbi_temp;
				for (j = 1; j < Pridictor_Length + 2; j++)
					Viterbi_error_Memory[i][j][0] = Viterbi_error_Memory[i / M_ary + Viterbi_temp * M_ary][j][1];
				if (y < Page_Size + 1 && y>0)	Viterbi_error_Memory[i][0][0] = eq_out[x][y - 1] - trellis[i][M_ary + Viterbi_temp];
			}

			// trace back and decide viterbi out
			if (y >= (Target_Order))
			{
				Viterbi_temp = 0;
				if (gg_pr_viterbi == 1)
				{
					for (i = 0; i < total_states; i++)
					{
						if ((i != 3) || (i != 13))
						{
							if (Viterbi_Path[Viterbi_temp][0] > Viterbi_Path[i][0])
								Viterbi_temp = i;
						}
					}
					if (y >= (depth - 1))
						viterbi_out[x][y - (depth - 1)] = int(Viterbi_State[Viterbi_temp][depth - 1][0]);

				}
				else
				{
					for (i = 0; i < total_states; i++)
						if (Viterbi_Path[Viterbi_temp][0] > Viterbi_Path[i][0])
							Viterbi_temp = i;
					if (y >= (depth - 1))
						viterbi_out[x][y - (depth - 1)] = int(Viterbi_State[Viterbi_temp][depth - 1][0]);
				}
			}

			// Shift the memory of path, state, and errors during the process
			for (i = 0; i < total_states; i++)
			{
				Viterbi_Path[i][1] = Viterbi_Path[i][0];
				for (j = depth - 1; j > 0; j--)
					Viterbi_State[i][j][1] = Viterbi_State[i][j - 1][0];
				for (j = depth - 1; j > 0; j--)
					Viterbi_whereState[i][j][1] = Viterbi_whereState[i][j - 1][0];
				for (j = 0; j < Pridictor_Length + 1; j++)
					Viterbi_error_Memory[i][j + 1][1] = Viterbi_error_Memory[i][j][0];
			}

			// prevent the path metric overflow
			// If all path metric is greater than 10000, all path metric is subtracted from the 10000.
			over = 0;
			for (i = 0; i < total_states; i++)
				if (Viterbi_Path[i][1] > 10000)	over++;

			if (over == total_states)
				for (i = 0; i < total_states; i++)	Viterbi_Path[i][1] -= 10000;
		}

		for (i = depth - 1; i >= 1; i--)
			viterbi_out[x][Page_Size - (i)] = int(Viterbi_State[Viterbi_temp][i][0]);

	}

	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			free(Viterbi_State[i][j]), free(Viterbi_whereState[i][j]);
	for (i = 0; i < total_states; i++)
		for (j = 0; j < Pridictor_Length + 2; j++)
			free(Viterbi_error_Memory[i][j]);
	for (i = 0; i < total_states; i++)
		free(Viterbi_State[i]), free(Viterbi_Path[i]), free(Viterbi_whereState[i]), free(Viterbi_error_Memory[i]);
	free(Viterbi_State), free(Viterbi_Path), free(Viterbi_whereState), free(Viterbi_error_Memory);
	free(differ);

}
void Hologram_SOVA_y(double** viterbi_out, double** eq_out, double** EI, double** trellis, double* Pridictor_Coef, int Pridictor_Length, int Target_Order, int Page_Size, int depth)
{
	Target_Order--;
	int i, j, k, x, y;				//// Temporary variables used during the for loop
	int Viterbi_temp;
	int vout = 0;
	double pri_error[M_ary];
	int over;
	double PM[M_ary];
	double alpha = 1.0;
	int error_viterbi = 0;
	int scale = 30;
	int iteration = 3;
	int total_states = (int)pow((double)M_ary, Target_Order);

	double eqout_total = 0, eqout_average;
	for (i = 0; i < Page_Size; i++)
		for (j = 0; j < Page_Size; j++)
			eqout_total += eq_out[i][j];
	eqout_average = eqout_total / (Page_Size * Page_Size);

	double** differ = (double**)malloc(total_states * sizeof(double*));
	for (i = 0; i < total_states; i++)
		differ[i] = (double*)calloc(M_ary, sizeof(double));

	for (i = 0; i < Page_Size; i++)
		for (j = 0; j < Page_Size; j++)
			viterbi_out[i][j] = 10000;

	unsigned long int*** Viterbi_State = (unsigned long int***)malloc(total_states * sizeof(unsigned long int**));
	for (i = 0; i < total_states; i++)
		Viterbi_State[i] = (unsigned long int**)malloc(depth * sizeof(unsigned long int*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			Viterbi_State[i][j] = (unsigned long int*)calloc(2, sizeof(unsigned long int));

	double** Viterbi_Path = (double**)malloc(total_states * sizeof(double*));
	for (i = 0; i < total_states; i++)
		Viterbi_Path[i] = (double*)calloc(2, sizeof(double));

	double*** Viterbi_whereState = (double***)malloc(total_states * sizeof(double**));
	for (i = 0; i < total_states; i++)
		Viterbi_whereState[i] = (double**)malloc(depth * sizeof(double*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			Viterbi_whereState[i][j] = (double*)calloc(2, sizeof(double));

	double*** Viterbi_error_Memory = (double***)malloc(total_states * sizeof(double**));
	for (i = 0; i < total_states; i++)
		Viterbi_error_Memory[i] = (double**)malloc((Pridictor_Length + 2) * sizeof(double*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < Pridictor_Length + 2; j++)
			Viterbi_error_Memory[i][j] = (double*)calloc(2, sizeof(double));

	for (y = 0; y < Page_Size; y++)
	{
		for (i = 1; i < total_states; i++)
			Viterbi_Path[i][1] = 10000, Viterbi_State[i][0][0] = 0, Viterbi_State[i][0][1] = 0, Viterbi_State[i][1][0] = 0, Viterbi_State[i][1][1] = 0;
		Viterbi_Path[0][1] = 0, Viterbi_State[0][0][0] = 0, Viterbi_State[0][0][1] = 0, Viterbi_State[0][1][0] = 0, Viterbi_State[0][1][1] = 0;

		for (i = 0; i < total_states; i++)
			for (j = 0; j < Pridictor_Length + 2; j++)
				Viterbi_error_Memory[i][j][0] = 0.0, Viterbi_error_Memory[i][j][1] = 0.0;

		for (x = 0; x < Page_Size + 1; x++)
		{
			for (i = 0; i < total_states; i++)
			{
				// calculate the predicted error
				pri_error[0] = 0.0, pri_error[1] = 0.0;
				for (j = 1; j < Pridictor_Length + 2; j++)
				{
					for (k = 0; k < M_ary; k++)
						pri_error[k] = Viterbi_error_Memory[i / M_ary + k * M_ary][j][1] * Pridictor_Coef[j - 1];
				}
				// calculate each path metric of a states
				if (NoiseFilter == 0)	for (k = 0; k < M_ary; k++) pri_error[k] = 0.0;
				if (x < Page_Size + 1 && x>0)
				{
					for (k = 0; k < M_ary; k++)
						PM[k] = Viterbi_Path[i / M_ary + k * M_ary][1] + pow(eq_out[x - 1][y] - trellis[i][M_ary + k] - pri_error[0], 2);// - (EI[x][y-1]/2)*(trellis[i][k]*2-1);
				}
				else
				{
					for (k = 0; k < M_ary; k++)
						PM[k] = Viterbi_Path[i / M_ary + k * M_ary][1];// +pow(eq_out[x][y - 1] - trellis[i][M_ary + k] - pri_error[0], 2);// - (EI[x][y-1]/2)*(trellis[i][k]*2-1);
				}
				// find minimum path of each state
				Viterbi_temp = 0;
				for (k = 0; k < M_ary; k++)
					if (PM[Viterbi_temp] > PM[k])	Viterbi_temp = k;

				// save differents between survivor path and competition path
				for (k = 0; k < M_ary; k++)
					differ[i][k] = PM[k] - PM[Viterbi_temp];

				// save the survivor path of each state
				Viterbi_Path[i][0] = PM[Viterbi_temp];
				for (j = 1; j < depth; j++)
					Viterbi_State[i][j][0] = Viterbi_State[i / M_ary + Viterbi_temp * M_ary][j][1];
				Viterbi_State[i][0][0] = (unsigned long)trellis[i][Viterbi_temp];
				for (j = 1; j < depth; j++)
					Viterbi_whereState[i][j][0] = Viterbi_whereState[i / M_ary + Viterbi_temp * M_ary][j][1];
				Viterbi_whereState[i][0][0] = Viterbi_temp;
				for (j = 1; j < Pridictor_Length + 2; j++)
					Viterbi_error_Memory[i][j][0] = Viterbi_error_Memory[i / M_ary + Viterbi_temp * M_ary][j][1];
				if (x < Page_Size + 1 && x>0)	Viterbi_error_Memory[i][0][0] = eq_out[x - 1][y] - trellis[i][M_ary + Viterbi_temp];
			}

			// trace back and decide viterbi out
			if (x >= (Target_Order))
			{
				Viterbi_temp = 0;
				for (i = 0; i < total_states; i++)
					if (Viterbi_Path[Viterbi_temp][0] > Viterbi_Path[i][0])	Viterbi_temp = i;

				if (x >= (depth - 1))
					viterbi_out[x - (depth - 1)][y] = int(Viterbi_State[Viterbi_temp][depth - 1][0]);
			}

			// Process for shifting the memory of path, state, and errors
			for (i = 0; i < total_states; i++)
			{
				Viterbi_Path[i][1] = Viterbi_Path[i][0];
				for (j = depth - 1; j > 0; j--)
					Viterbi_State[i][j][1] = Viterbi_State[i][j - 1][0];
				for (j = depth - 1; j > 0; j--)
					Viterbi_whereState[i][j][1] = Viterbi_whereState[i][j - 1][0];
				for (j = 0; j < Pridictor_Length + 1; j++)
					Viterbi_error_Memory[i][j + 1][1] = Viterbi_error_Memory[i][j][0];
			}

			// prevent the path metric overflow
			// If all path metric is greater than 10000, all path metric is subtracted from the 10000.
			over = 0;
			for (i = 0; i < total_states; i++)
				if (Viterbi_Path[i][1] > 10000)	over++;

			if (over == total_states)
				for (i = 0; i < total_states; i++)	Viterbi_Path[i][1] -= 10000;
		}

		for (i = depth - 1; i >= 1; i--)
			viterbi_out[Page_Size - (i)][y] = int(Viterbi_State[Viterbi_temp][i][0]);
	}

	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			free(Viterbi_State[i][j]), free(Viterbi_whereState[i][j]);
	for (i = 0; i < total_states; i++)
		for (j = 0; j < Pridictor_Length + 2; j++)
			free(Viterbi_error_Memory[i][j]);
	for (i = 0; i < total_states; i++)
		free(Viterbi_State[i]), free(Viterbi_Path[i]), free(Viterbi_whereState[i]), free(Viterbi_error_Memory[i]);
	free(Viterbi_State), free(Viterbi_Path), free(Viterbi_whereState), free(Viterbi_error_Memory);
	free(differ);
}
// 7 To combat ISI, based on the one-dimensional (1D) partial response maximum likelihood (PRML) detection scheme
void PRML_Convolution_gg_grad_23(double** trellis, double* Target_Coef, int Target_Order, int in_size, int* input)
{
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//// trellis[i][j]	i=state, j= 0,1 => reference_bit(0 up, 1 down) 2,3 => reference_value(2 up, 3 down) ////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Target_Order--;
	int i, j, k;					//// Temporary variables used during the for loop
	int total_states = (int)pow((double)in_size, Target_Order);
	int ggg = 0;
	int* temp_trel_a = (int*)calloc(Target_Order, sizeof(int));


	//// setting reference_bit ////
	for (i = 0; i < total_states; i++)
	{
		for (j = 0; j < in_size; j++)
			trellis[i][j] = input[i % in_size];
	}
	///////////////////////////////

	for (i = 0; i < Target_Order; i++)
		temp_trel_a[i] = input[0];										// Initialize to 0
	for (i = 0; i < total_states; i++)
	{
		temp_trel_a[0] = input[i % in_size];
		temp_trel_a[1] = input[i / in_size];

		for (k = 0; k < in_size; k++)
		{
			trellis[i][in_size + k] = 0;
			for (j = 0; j < Target_Order; j++)
				trellis[i][in_size + k] += Target_Coef[j] * temp_trel_a[j];
			trellis[i][in_size + k] += (Target_Coef[Target_Order] * input[k]);
			trellis[i][in_size + k] /= M_ary - 1;

		}
	}

	free(temp_trel_a);

}
void PRML_Convolution_for_modul_gg_grad_23(int*** trellis_modul)
{
	int i, j, k;
	int total_states = (int)pow((double)M_ary, 2);
	int temp_piece;
	int ggg = 0;

	int odd_piece[16][3] = { {0,1,1},{0,1,2},{0,1,3},{0,2,1},{0,2,2},{0,2,3},{1,1,1},{1,1,2},
						  {1,1,3},{1,2,1},{1,2,2},{1,2,3},{2,1,1},{2,1,2},{2,1,3},{2,2,1} };

	for (i = 0; i < total_states; i++)
	{
		for (j = 0; j < total_states; j++)
		{
			temp_piece = (i + j) % 16;
			for (k = 0; k < 3; k++)
				trellis_modul[i][j][k] = odd_piece[temp_piece][k];
			ggg = 0;
		}
	}
}
void PRML_Convolution_for_modul_gg_grad_23_prof(int*** trellis_modul)
{

	int i, j, k;
	int total_states = (int)pow((double)M_ary, 2);
	int temp_piece;
	int ggg = 0;

	int piece[16][3] = { {0,0,1},{0,0,2},{0,0,3},{0,1,1},{0,1,2},{0,1,3},{0,2,2},{0,2,3},
					  {0,3,3},{1,1,2},{1,1,3},{1,2,2},{1,2,3},{1,3,3},{2,2,3},{2,3,3} };

	for (i = 0; i < total_states; i++)
	{
		for (j = 0; j < total_states; j++)
		{
			temp_piece = (i + j) % 16;
			for (k = 0; k < 3; k++)
				trellis_modul[i][j][k] = piece[temp_piece][k];
			ggg = 0;
		}
	}
}
void PRML_Convolution_for_modul_gg_dmin2(int*** trellis_modul)
{

	int i, j, k;
	int total_states = (int)pow((double)M_ary, 2);
	int temp_piece;
	int ggg = 0;

	int piece[16][3] = { {0,0,1},{1,3,3},{0,0,3},{1,3,1},{0,1,0},{1,2,2},{0,1,2},{1,2,0},{0,2,1},{1,1,3},{0,2,3},{1,1,1},{0,3,0},{1,0,2},{0,3,2},{1,0,0} };

	for (i = 0; i < total_states; i++)
	{
		for (j = 0; j < total_states; j++)
		{
			temp_piece = (i + j) % 16;
			for (k = 0; k < 3; k++)
				trellis_modul[i][j][k] = piece[temp_piece][k];
			ggg = 0;
		}
	}
}
void PRML_Convolution_for_modul_gg_grad_23_prof_method0(int**** trellis_modul)
{

	int i, j, k, l;
	int total_states = ((int)pow((double)M_ary, 2)) / 2;
	int temp_piece;
	int ggg = 0;

	int piece[16][3] = { {0,0,1},{0,0,2},{0,0,3},{0,1,1},{0,1,2},{0,1,3},{0,2,2},{0,2,3},
					  {0,3,3},{1,1,2},{1,1,3},{1,2,2},{1,2,3},{1,3,3},{2,2,3},{2,3,3} };


	for (i = 0; i < total_states; i++)
	{
		for (j = 0; j < total_states; j++)
		{
			for (k = 0; k < 2; k++)
			{
				temp_piece = ((i * 2) + (j * 2) + k) % 16;
				ggg = 0;
				for (l = 0; l < 3; l++)
					trellis_modul[i][j][k][l] = piece[temp_piece][l];
			}
		}
	}
}
void PRML_Convolution_for_modul_gg_dmin2_method0(int**** trellis_modul)
{

	int i, j, k, l;
	int total_states = ((int)pow((double)M_ary, 2)) / 2;
	int temp_piece;
	int ggg = 0;

	int piece[16][3] = { {0,0,1},{1,3,3},{0,0,3},{1,3,1},{0,1,0},{1,2,2},{0,1,2},{1,2,0},{0,2,1},{1,1,3},{0,2,3},{1,1,1},{0,3,0},{1,0,2},{0,3,2},{1,0,0} };

	for (i = 0; i < total_states; i++)
	{
		for (j = 0; j < total_states; j++)
		{
			for (k = 0; k < 2; k++)
			{
				temp_piece = ((i * 2) + (j * 2) + k) % 16;
				ggg = 0;
				for (l = 0; l < 3; l++)
					trellis_modul[i][j][k][l] = piece[temp_piece][l];
			}
		}
	}
}

void Hologram_SOVA_x_for_gg_grad(double** viterbi_out, double** eq_out, double** EI, double** trellis, double** trellisA, double* Pridictor_Coef, int Pridictor_Length, int Target_Order, int Page_Size, int depth)
{
	Target_Order--;
	int i, j, k, x, y;				//// Temporary variables used during the for loop
	int Viterbi_temp;
	int vout = 0;
	double pri_error[M_ary];
	int over;
	double PM[M_ary];
	double alpha = 1.0;
	int error_viterbi = 0;
	int scale = 30;
	int iteration = 3;
	int total_states = (int)pow((double)M_ary, Target_Order);
	int test_line = 0;
	int ggg = 0;


	double eqout_total = 0, eqout_average;
	for (i = 0; i < Page_Size; i++)
		for (j = 0; j < Page_Size; j++)
			eqout_total += eq_out[i][j];
	eqout_average = eqout_total / (Page_Size * Page_Size);

	double** differ = (double**)malloc(total_states * sizeof(double*));
	for (i = 0; i < total_states; i++)
		differ[i] = (double*)calloc(M_ary, sizeof(double));

	for (i = 0; i < Page_Size; i++)
		for (j = 0; j < Page_Size; j++)
			viterbi_out[i][j] = 10000;

	unsigned long int*** Viterbi_State = (unsigned long int***)malloc(total_states * sizeof(unsigned long int**));
	for (i = 0; i < total_states; i++)
		Viterbi_State[i] = (unsigned long int**)malloc(depth * sizeof(unsigned long int*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			Viterbi_State[i][j] = (unsigned long int*)calloc(2, sizeof(unsigned long int));

	double** Viterbi_Path = (double**)malloc(total_states * sizeof(double*));
	for (i = 0; i < total_states; i++)
		Viterbi_Path[i] = (double*)calloc(2, sizeof(double));

	double*** Viterbi_whereState = (double***)malloc(total_states * sizeof(double**));
	for (i = 0; i < total_states; i++)
		Viterbi_whereState[i] = (double**)malloc(depth * sizeof(double*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			Viterbi_whereState[i][j] = (double*)calloc(2, sizeof(double));

	double*** Viterbi_error_Memory = (double***)malloc(total_states * sizeof(double**));
	for (i = 0; i < total_states; i++)
		Viterbi_error_Memory[i] = (double**)malloc((Pridictor_Length + 2) * sizeof(double*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < Pridictor_Length + 2; j++)
			Viterbi_error_Memory[i][j] = (double*)calloc(2, sizeof(double));

	for (x = 0; x < Page_Size; x++)
	{
		for (i = 1; i < total_states; i++)
			Viterbi_Path[i][1] = 10000, Viterbi_State[i][0][0] = 0, Viterbi_State[i][0][1] = 0, Viterbi_State[i][1][0] = 0, Viterbi_State[i][1][1] = 0;
		Viterbi_Path[0][1] = 0, Viterbi_State[0][0][0] = 0, Viterbi_State[0][0][1] = 0, Viterbi_State[0][1][0] = 0, Viterbi_State[0][1][1] = 0;

		for (i = 0; i < total_states; i++)
			for (j = 0; j < Pridictor_Length + 2; j++)
				Viterbi_error_Memory[i][j][0] = 0.0, Viterbi_error_Memory[i][j][1] = 0.0;


		test_line = x % 3;
		if (test_line == 0)
		{
			total_states = (int)pow((double)2, Target_Order);

			for (y = 0; y < Page_Size + 1; y++)
			{
				ggg = 0;
				for (i = 0; i < total_states; i++)
				{
					// calculate the predicted error
					pri_error[0] = 0.0, pri_error[1] = 0.0;
					for (j = 1; j < Pridictor_Length + 2; j++)
					{
						for (k = 0; k < 2; k++)
							pri_error[k] = Viterbi_error_Memory[i / 2 + k * 2][j][1] * Pridictor_Coef[j - 1];
					}
					// calculate each path metric of a states
					if (NoiseFilter == 0)	for (k = 0; k < 2; k++) pri_error[k] = 0.0;

					if (y < Page_Size + 1 && y>0)
					{
						for (k = 0; k < 2; k++)
							PM[k] = Viterbi_Path[i / 2 + k * 2][1] + pow(eq_out[x][y - 1] - trellisA[i][2 + k] - pri_error[0], 2);// - (EI[x][y-1]/2)*(trellis[i][k]*2-1);
					}
					else
					{
						for (k = 0; k < 2; k++)
							PM[k] = Viterbi_Path[i / 2 + k * 2][1];//+pow(eq_out[x][y-1]-trellis[i][M_ary+k]-pri_error[0],2);// - (EI[x][y-1]/2)*(trellis[i][k]*2-1);
					}

					// find minimum path of each state				
					Viterbi_temp = 0;
					for (k = 0; k < 2; k++)
						if (PM[Viterbi_temp] > PM[k])
							Viterbi_temp = k;

					// save differents between survivor path and competition path
					for (k = 0; k < 2; k++)
						differ[i][k] = PM[k] - PM[Viterbi_temp];

					// save the survivor path of each state
					Viterbi_Path[i][0] = PM[Viterbi_temp];
					for (j = 1; j < depth; j++)
						Viterbi_State[i][j][0] = Viterbi_State[i / 2 + Viterbi_temp * 2][j][1];
					Viterbi_State[i][0][0] = (unsigned long)trellisA[i][Viterbi_temp];
					for (j = 1; j < depth; j++)
						Viterbi_whereState[i][j][0] = Viterbi_whereState[i / 2 + Viterbi_temp * 2][j][1];
					Viterbi_whereState[i][0][0] = Viterbi_temp;
					for (j = 1; j < Pridictor_Length + 2; j++)
						Viterbi_error_Memory[i][j][0] = Viterbi_error_Memory[i / 2 + Viterbi_temp * 2][j][1];
					if (y < Page_Size + 1 && y>0)	Viterbi_error_Memory[i][0][0] = eq_out[x][y - 1] - trellisA[i][2 + Viterbi_temp];
				}

				// trace back and decide viterbi out
				if (y >= (Target_Order))
				{
					Viterbi_temp = 0;
					for (i = 0; i < total_states; i++)
						if (Viterbi_Path[Viterbi_temp][0] > Viterbi_Path[i][0])
							Viterbi_temp = i;
					if (y >= (depth - 1))
						viterbi_out[x][y - (depth - 1)] = int(Viterbi_State[Viterbi_temp][depth - 1][0]);
				}

				// shift the memory of path and state and errors¸Þ¸ð¸®¸¦ shift½ÃÅ°´Â ±¸Á¶
				for (i = 0; i < total_states; i++)
				{
					Viterbi_Path[i][1] = Viterbi_Path[i][0];
					for (j = depth - 1; j > 0; j--)
						Viterbi_State[i][j][1] = Viterbi_State[i][j - 1][0];
					for (j = depth - 1; j > 0; j--)
						Viterbi_whereState[i][j][1] = Viterbi_whereState[i][j - 1][0];
					for (j = 0; j < Pridictor_Length + 1; j++)
						Viterbi_error_Memory[i][j + 1][1] = Viterbi_error_Memory[i][j][0];
				}

				// prevent the path metric overflow
				// If all path metric is greater than 10000, all path metric is subtracted from the 10000.
				over = 0;
				for (i = 0; i < total_states; i++)
					if (Viterbi_Path[i][1] > 10000)	over++;

				if (over == total_states)
					for (i = 0; i < total_states; i++)	Viterbi_Path[i][1] -= 10000;
			}
		}
		else
		{
			total_states = (int)pow((double)4, Target_Order);

			for (y = 0; y < Page_Size + 1; y++)
			{
				ggg = 0;
				for (i = 0; i < total_states; i++)
				{
					// calculate the predicted error
					pri_error[0] = 0.0, pri_error[1] = 0.0;
					for (j = 1; j < Pridictor_Length + 2; j++)
					{
						for (k = 0; k < M_ary; k++)
							pri_error[k] = Viterbi_error_Memory[i / M_ary + k * M_ary][j][1] * Pridictor_Coef[j - 1];
					}
					// calculate each path metric of a states
					if (NoiseFilter == 0)	for (k = 0; k < M_ary; k++) pri_error[k] = 0.0;

					if (y < Page_Size + 1 && y>0)
					{
						for (k = 0; k < M_ary; k++)
							PM[k] = Viterbi_Path[i / M_ary + k * M_ary][1] + pow(eq_out[x][y - 1] - trellis[i][M_ary + k] - pri_error[0], 2);// - (EI[x][y-1]/2)*(trellis[i][k]*2-1);
					}
					else
					{
						for (k = 0; k < M_ary; k++)
							PM[k] = Viterbi_Path[i / M_ary + k * M_ary][1];//+pow(eq_out[x][y-1]-trellis[i][M_ary+k]-pri_error[0],2);// - (EI[x][y-1]/2)*(trellis[i][k]*2-1);
					}

					// find minimum path of each state		

					Viterbi_temp = 0;
					for (k = 0; k < M_ary; k++)
						if (PM[Viterbi_temp] > PM[k])
							Viterbi_temp = k;

					// save differents between survivor path and competition path
					for (k = 0; k < M_ary; k++)
						differ[i][k] = PM[k] - PM[Viterbi_temp];

					// save the survivor path of each state
					Viterbi_Path[i][0] = PM[Viterbi_temp];
					for (j = 1; j < depth; j++)
						Viterbi_State[i][j][0] = Viterbi_State[i / M_ary + Viterbi_temp * M_ary][j][1];
					Viterbi_State[i][0][0] = (unsigned long)trellis[i][Viterbi_temp];
					for (j = 1; j < depth; j++)
						Viterbi_whereState[i][j][0] = Viterbi_whereState[i / M_ary + Viterbi_temp * M_ary][j][1];
					Viterbi_whereState[i][0][0] = Viterbi_temp;
					for (j = 1; j < Pridictor_Length + 2; j++)
						Viterbi_error_Memory[i][j][0] = Viterbi_error_Memory[i / M_ary + Viterbi_temp * M_ary][j][1];
					if (y < Page_Size + 1 && y>0)	Viterbi_error_Memory[i][0][0] = eq_out[x][y - 1] - trellis[i][M_ary + Viterbi_temp];
				}

				// trace back and decide viterbi out
				if (y >= (Target_Order))
				{
					Viterbi_temp = 0;
					for (i = 0; i < total_states; i++)
						if (Viterbi_Path[Viterbi_temp][0] > Viterbi_Path[i][0])
							Viterbi_temp = i;
					if (y >= (depth - 1))
						viterbi_out[x][y - (depth - 1)] = int(Viterbi_State[Viterbi_temp][depth - 1][0]);
				}


				// Logic for shifting the memory of path, state, and errors
				for (i = 0; i < total_states; i++)
				{
					Viterbi_Path[i][1] = Viterbi_Path[i][0];
					for (j = depth - 1; j > 0; j--)
						Viterbi_State[i][j][1] = Viterbi_State[i][j - 1][0];
					for (j = depth - 1; j > 0; j--)
						Viterbi_whereState[i][j][1] = Viterbi_whereState[i][j - 1][0];
					for (j = 0; j < Pridictor_Length + 1; j++)
						Viterbi_error_Memory[i][j + 1][1] = Viterbi_error_Memory[i][j][0];
				}

				// prevent the path metric overflow
				// If all path metric is greater than 10000, all path metric is subtracted from the 10000.
				over = 0;
				for (i = 0; i < total_states; i++)
					if (Viterbi_Path[i][1] > 10000)	over++;


				if (over == total_states)
					for (i = 0; i < total_states; i++)	Viterbi_Path[i][1] -= 10000;
			}

		}

		for (i = depth - 1; i >= 1; i--)
			viterbi_out[x][Page_Size - (i)] = int(Viterbi_State[Viterbi_temp][i][0]);


	}
	total_states = (int)pow((double)4, Target_Order);

	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			free(Viterbi_State[i][j]), free(Viterbi_whereState[i][j]);
	for (i = 0; i < total_states; i++)
		for (j = 0; j < Pridictor_Length + 2; j++)
			free(Viterbi_error_Memory[i][j]);
	for (i = 0; i < total_states; i++)
		free(Viterbi_State[i]), free(Viterbi_Path[i]), free(Viterbi_whereState[i]), free(Viterbi_error_Memory[i]);
	free(Viterbi_State), free(Viterbi_Path), free(Viterbi_whereState), free(Viterbi_error_Memory);
	free(differ);

}

