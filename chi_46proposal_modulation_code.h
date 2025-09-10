#include <iostream>
void Encode_chiProposal46_4ary(int** PAGE, int* input_1Ddata, int Page_Size)
{
	////////////////////////////////////////////////////////////////////////
	//		simple encoding scheme 4 symbol (8 bit) -> 6 symbol (12 bit)  //
	//																	  //
	////////////////////////////////////////////////////////////////////////

	int i, j, k;
	int input[4];
	int temp_piece;

	// Bao nhiêu phần từ ở hàng dọc => i là index hàng dọc
	for (i = 0; i < Page_Size / 3; i++)
	{
		// Bao nhiêu phần tử ở hàng ngang => j là index hàng ngang
		for (j = 0; j < Page_Size / 2; j++)
		{
			//temp_piece = 0;
			for (k = 0; k < 4; k++)
			{
				// Chuyển hóa từng bit 1 trong 4 bit đầu vào
				input[k] = input_1Ddata[(i * Page_Size / 2 + j) * 4 + k];
				//temp_piece += input[k] << (3 - k);
			}
			// x0
			PAGE[i * 3 + 0][j * 2 + 0] = input[0];
			// x5
			PAGE[i * 3 + 2][j * 2 + 1] = input[1];

			// x1 = 1 + mod(a/2,2)
			PAGE[i * 3 + 0][j * 2 + 1] = (input[2] >> 1) % 2 + 1;
			// x2
			PAGE[i * 3 + 1][j * 2 + 0] = input[2] % 2 + 1;

			// x3
			PAGE[i * 3 + 1][j * 2 + 1] = (input[3] >> 1) % 2 + 1;
			// x4
			PAGE[i * 3 + 2][j * 2 + 0] = input[3] % 2 + 1;

			//temp_piece = input[0] << 2;
			//temp_piece += input[1] << 0;			
			//for (k = 0; k < 4; k++)//so bit dau ra
			//	PAGE[i * 4 + k][j] = piece45[temp_piece][k];
		}
	}
}
void Decode_chiProposal46_4ary(int* output_1Ddata, double** PAGE, int Page_Size)
{
	//std::cout << "10x10 before decode:" << std::endl;
	//for (int i = 0; i < 10 && i < Page_Size; i++) {
	//	printf("Row %2d: ", i + 1);
	//	for (int j = 0; j < 10 && j < Page_Size; j++) { // Chỉ log tối đa 10 cột
	//		printf("%3d ", (int)PAGE[i][j]);
	//	}
	//	printf("\n");
	//}

	//for (int i = 0; i < 10 && i < Page_Size; i++) {
	//	printf("Row %2d: ", i + 1);
	//	for (int j = 0; j < 10 && j < Page_Size; j++) { // Chỉ log tối đa 10 cột
	//		printf("%3.2f ", PAGE[i][j]);
	//	}
	//	printf("\n");
	//}
	//////////////////////////////////////////////////////////////////////
	//		simple encoding scheme 3 symbol (6bit) -> 4 symbol (8bit)	//
	//	
	//////////////////////////////////////////////////////////////////////

	int i, j, k, m, l;
	int temp;

	for (i = 0; i < Page_Size / 3; i++)
	{
		for (j = 0; j < Page_Size / 2; j++)
		{
			output_1Ddata[(i * Page_Size / 2 + j) * 4 + 0] = int(PAGE[i * 3 + 0][j * 2 + 0]);
			output_1Ddata[(i * Page_Size / 2 + j) * 4 + 1] = int(PAGE[i * 3 + 2][j * 2 + 1]);

			output_1Ddata[(i * Page_Size / 2 + j) * 4 + 2] = int(PAGE[i * 3 + 0][j * 2 + 1] > 1 ? 1 : 0) * 2 + int(PAGE[i * 3 + 1][j * 2 + 0] > 1 ? 1 : 0);
			output_1Ddata[(i * Page_Size / 2 + j) * 4 + 3] = int(PAGE[i * 3 + 1][j * 2 + 1] > 1 ? 1 : 0) * 2 + int(PAGE[i * 3 + 2][j * 2 + 0] > 1 ? 1 : 0);
		}
	}

	/*
	double BM[16] = { 0.0, };

	for (i = 0; i < Page_Size / 4; i++)
	{
		for (j = 0; j < Page_Size; j++)
		{
			int b[2] = { 0, };
			temp = 0;
			for (l = 0; l < 16; l++)
			{
				BM[l] = 0;
				for (k = 0; k < 4; k++)//so bit dau ra
					BM[l] += pow((PAGE[i * 4 + k][j] - piece45[l][k]), 2);

				if (BM[temp] > BM[l])
					temp = l;
			}

			//deci2Mary(temp, 4, b);
			//for (m = 1; m >= 0; m--)
			//{
			//	output_1Ddata[(i*Page_Size + j) * 2 + m] = temp % 4;// , printf("%d ", output_1Ddata[(i*Page_Size + j) * 4 + m]);
			//	temp = temp / 4;
			//}
			output_1Ddata[(i*Page_Size + j) * 2 + 0] = temp / 4;
			output_1Ddata[(i*Page_Size + j) * 2 + 1] = temp % 4;

			//printf("\n");
			//for (m = 0; m < 4; m++)
				//output_1Ddata[(i*Page_Size + j) * 4 + m] = b[m];
		}
	}
	*/
}
