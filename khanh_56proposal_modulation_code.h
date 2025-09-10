#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <iomanip> 
void printMatrix(int** matrix, int height, int width);

void encodeBits(int* unencodedBits, int unencodedSize, int*& encodedBits, int& encodedSize) {
	int groupSize = 5; // Each group of 5 bits will be encoded into 6 bits
	int encodedGroupSize = 6;
	int numGroups = unencodedSize / groupSize;

	encodedSize = numGroups * encodedGroupSize;
	encodedBits = new int[encodedSize];
	//std::cout << "Encoding process:\n";
	for (int i = 0; i < numGroups; ++i) {
		//std::cout << "Group " << i + 1 << " (original bits): ";
		// Copy the first 5 bits
		for (int j = 0; j < groupSize; ++j) {
			encodedBits[i * encodedGroupSize + j] = unencodedBits[i * groupSize + j];
			//std::cout << unencodedBits[i * groupSize + j] << " ";
		}
		// Add a random 6th bit (values from 0 to 3)
		encodedBits[i * encodedGroupSize + groupSize] = rand() % 4;

		//// print encoded 
		//std::cout << " | Encoded group: ";
		//for (int j = 0; j < encodedGroupSize; ++j) {
		//	std::cout << encodedBits[i * encodedGroupSize + j] << " ";
		//}
		//std::cout << "\n";
	}
}

void mapBitsToMatrix(int** matrix, int* bits, int matrixRows, int matrixCols) {
	// Divide the bits array into two parts
	int aBits[6] = { bits[0], bits[1], bits[2], bits[3], bits[4], bits[5] }; // aBits
	int xyBits[6] = { bits[11], bits[10], bits[9], bits[8], bits[7], bits[6] }; // xyBits (reversed order)

	// Place xyBits into the first row of the matrix
	int xyIndex = 0;
	for (int j = 0; j < matrixCols; ++j) {
		matrix[0][j] = xyBits[xyIndex++];
	}

	// Place aBits and the remaining xyBits into subsequent rows
	int aIndex = 0;
	for (int i = 1; i < matrixRows; ++i) {
		// Fill in aBits first
		for (int j = 0; j < i; ++j) { // Fill a number of bits corresponding to the row index
			if (aIndex < 6) { // Ensure we don't exceed the size of aBits
				matrix[i][j] = aBits[aIndex++];
			}
		}

		// Fill in the remaining xyBits
		for (int j = i; j < matrixCols; ++j) {
			if (xyIndex < 6) { // Ensure we don't exceed the size of xyBits
				matrix[i][j] = xyBits[xyIndex++];
			}
		}
	}
	//printMatrix(matrix, matrixRows, matrixCols);
}

void printMatrix(int** matrix, int height, int width) {
	if (!matrix) {
		std::cerr << "Error: Matrix is null!" << std::endl;
		return;
	}

	int cellWidth = 4; // Độ rộng của mỗi ô
	for (int i = 0; i < height; ++i) {
		if (!matrix[i]) { // Kiểm tra hàng chưa được cấp phát
			std::cerr << "Error: Row " << i << " is null!" << std::endl;
			continue;
		}

		for (int j = 0; j < width; ++j) {
			std::cout << std::setw(cellWidth) << matrix[i][j]; // Sử dụng setw để căn lề
		}
		std::cout << std::endl;
	}
}

void Encode_khanhProposal_5_6_4ary(int** PAGE, int* unencodedBits, int Page_Size) {
	int Page_Height = Page_Size;
	int Page_Width = Page_Size;
	int Block_Height = 4, Block_Width = 3, Bits_Per_Block = Block_Height * Block_Width;

	int* encodedBits = nullptr;
	int encodedSize = 0;
	encodeBits(unencodedBits, Page_Size, encodedBits, encodedSize);

	//// Log 100 phần tử đầu tiên của unencodedBits
	//printf("=== Logging first 50 unencodedBits ===\n");
	//for (int i = 0; i < 50; i++) {
	//	if (i % 5 == 0 && i != 0) printf("\n"); // Xuống dòng sau mỗi 10 phần tử để dễ nhìn
	//	printf("%d ", unencodedBits[i]);
	//}
	//printf("\n\n");

	// Log 100 phần tử đầu tiên của encodedBits
	//printf("=== Logging first 60 encodedBits ===\n");
	//for (int i = 0; i < 30 && i < encodedSize; i++) { // Đảm bảo không vượt quá kích thước mảng
	//	if (i % 6 == 0 && i != 0) printf("\n"); // Xuống dòng sau mỗi 10 phần tử
	//	printf("%d ", encodedBits[i]);
	//}
	//printf("\n");

	//return;


	int bitIndex = 0, currentRow = 0, currentCol = 0;

	// Cấp phát bộ nhớ cho block chỉ một lần
	int** block = new int* [Block_Height];
	for (int i = 0; i < Block_Height; ++i) {
		block[i] = new int[Block_Width];
	}

	int* currentBits = new int[Bits_Per_Block]; // Dùng chung cho tất cả vòng lặp

	while (bitIndex < encodedSize) {
		// Reset lại giá trị block thay vì cấp phát lại
		for (int i = 0; i < Block_Height; ++i) {
			std::fill(block[i], block[i] + Block_Width, 0); // Sử dụng fill để đặt 0 nhanh hơn
		}
		/*std::cout << bitIndex << std::endl;
		std::cout << encodedSize << std::endl;*/
		if (bitIndex + Bits_Per_Block > encodedSize) break;
		if (bitIndex + Bits_Per_Block <= encodedSize) {
			// Copy bits vào mảng dùng chung
			memcpy(currentBits, encodedBits + bitIndex, Bits_Per_Block * sizeof(int));
			bitIndex += Bits_Per_Block;

			// Ánh xạ bits vào ma trận block
			mapBitsToMatrix(block, currentBits, Block_Height, Block_Width);

			// Tính offsetCol
			int offsetCol = ((currentRow / Block_Height) % 2 == 1) ? 1 : 0;

			for (int i = 0; i < Block_Height; ++i) {
				for (int j = 0; j < Block_Width; ++j) {
					int targetCol = currentCol + j + offsetCol;
					if (offsetCol > 0 && currentCol == 0 && j == 0) {
						PAGE[currentRow + i][currentCol + j] = 1; // Gán giá trị mặc định
					}
					if (targetCol < Page_Width) {
						PAGE[currentRow + i][targetCol] = block[i][j];
					}
				}
			}
			if (offsetCol > 0 && currentRow > 0) {
				for (int i = 0; i < Block_Height; i++) {
					PAGE[currentRow - Block_Height + i][Page_Width - 1] = 1; // Điền 1 vào cột cuối của hàng trước
				}
			}
			currentCol += Block_Width;

			if (currentCol + Block_Width > Page_Width) {
				currentRow += Block_Height;
				currentCol = 0;
			}

			/*std::cout << "Done" << std::endl;*/
		}
	}
	//std::cout << "10x10 after encoded:" << std::endl;
	//for (int i = 0; i < 10 && i < Page_Size; i++) {
	//	printf("Row %2d: ", i + 1);
	//	for (int j = 0; j < 10 && j < Page_Size; j++) { // Chỉ log tối đa 10 cột
	//		printf("%3d ", (int)PAGE[i][j]);
	//	}
	//	printf("\n");
	//}

	// Giải phóng bộ nhớ block và currentBits sau vòng lặp
	for (int i = 0; i < Block_Height; ++i) {
		delete[] block[i];
	}
	delete[] block;
	delete[] currentBits;

	delete[] encodedBits;
}




void Decode_khanhProposal_5_6_4ary(int* decodedBitsArray, double** page, int Page_Size) {
	//std::cout << "10x10 before decode:" << std::endl;
	//for (int i = 0; i < 10 && i < Page_Size; i++) {
	//	printf("Row %2d: ", i + 1);
	//	for (int j = 0; j < 10 && j < Page_Size; j++) { // Chỉ log tối đa 10 cột
	//		printf("%3d ", (int)page[i][j]);
	//	}
	//	printf("\n");
	//}

	int Page_Height = Page_Size;
	int Page_Width = Page_Size;
	int Block_Height = 4; // Height of each block
	int Block_Width = 3;  // Width of each block
	int Bits_Per_Block = Block_Height * Block_Width; // Total bits per block

	std::vector<int> decodedBits; // Temporary container for decoded bits

	for (int currentRow = 0; currentRow < Page_Height; currentRow += Block_Height) {
		int offsetCol = ((currentRow / Block_Height) % 2 == 1) ? 1 : 0;

		for (int currentCol = 0; currentCol < Page_Width; currentCol += Block_Width) {
			// Ensure the block is within the bounds of the page
			if (currentCol + offsetCol + Block_Width - 1 < Page_Width && currentRow + Block_Height - 1 < Page_Height) {
				std::vector<int> blockBits(Bits_Per_Block - 2, 0);

				// Extract bits from the page into blockBits (convert double to int explicitly)
				blockBits[0] = static_cast<int>(page[currentRow + 1][currentCol + offsetCol]);          // (1,0)
				blockBits[1] = static_cast<int>(page[currentRow + 2][currentCol + offsetCol]);          // (2,0)
				blockBits[2] = static_cast<int>(page[currentRow + 2][currentCol + 1 + offsetCol]);      // (2,1)
				blockBits[3] = static_cast<int>(page[currentRow + 3][currentCol + offsetCol]);          // (3,0)
				blockBits[4] = static_cast<int>(page[currentRow + 3][currentCol + 1 + offsetCol]);      // (3,1)

				// Extract xyBits and reverse them to match the original logic
				std::vector<int> xyBits = {
					static_cast<int>(page[currentRow][currentCol + 1 + offsetCol]),                     // (0,1)
					static_cast<int>(page[currentRow][currentCol + 2 + offsetCol]),                     // (0,2)
					static_cast<int>(page[currentRow + 1][currentCol + 1 + offsetCol]),                 // (1,1)
					static_cast<int>(page[currentRow + 1][currentCol + 2 + offsetCol]),                 // (1,2)
					static_cast<int>(page[currentRow + 2][currentCol + 2 + offsetCol])                  // (2,2)
				};

				std::reverse(xyBits.begin(), xyBits.end());

				// Add reversed xyBits to blockBits
				for (int i = 0; i < 5; ++i) {
					blockBits[5 + i] = xyBits[i];
				}

				// Append blockBits to decodedBits
				decodedBits.insert(decodedBits.end(), blockBits.begin(), blockBits.end());
			}
		}
	}

	// Copy decodedBits back into decodedBitsArray
	for (size_t i = 0; i < decodedBits.size(); ++i) {
		decodedBitsArray[i] = decodedBits[i];
	}
}


void Encode_khanhProposal_5_6_4ary1(int** PAGE, int* unencodedBits, int Page_Size) {
	int Page_Height = Page_Size;
	int Page_Width = Page_Size;
	int Block_Height = 4, Block_Width = 3, Bits_Per_Block = Block_Height * Block_Width;

	int* encodedBits = nullptr;
	int encodedSize = 0;
	encodeBits(unencodedBits, Page_Size, encodedBits, encodedSize);

	int bitIndex = 0, currentRow = 0, currentCol = 0;

	while (bitIndex < encodedSize) {
		int** block = new int* [Block_Height];
		for (int i = 0; i < Block_Height; ++i) {
			block[i] = new int[Block_Width];
			memset(block[i], 0, Block_Width * sizeof(int));
		}
		if (bitIndex + Bits_Per_Block > encodedSize) break;
		if (bitIndex + Bits_Per_Block <= encodedSize) {
			int* currentBits = new int[Bits_Per_Block];
			memcpy(currentBits, encodedBits + bitIndex, Bits_Per_Block * sizeof(int));
			bitIndex += Bits_Per_Block;

			mapBitsToMatrix(block, currentBits, Block_Height, Block_Width);

			int offsetCol = ((currentRow / Block_Height) % 2 == 1) ? 1 : 0;

			for (int i = 0; i < Block_Height; ++i) {
				for (int j = 0; j < Block_Width; ++j) {
					int targetCol = currentCol + j + offsetCol;
					if (offsetCol > 0 && currentCol == 0 && j == 0) {
						PAGE[currentRow + i][currentCol + j] = 1; // Gán giá trị mặc định
					}
					if (targetCol < Page_Width) {
						PAGE[currentRow + i][targetCol] = block[i][j];
					}
				}
			}
			if (offsetCol > 0 && currentRow > 0) {
				for (int i = 0; i < Block_Height; i++) {
					PAGE[currentRow - Block_Height + i][Page_Width - 1] = 1; // Điền 1 vào cột cuối của hàng trước
				}
			}
			currentCol += Block_Width;

			if (currentCol + Block_Width > Page_Width) {
				currentRow += Block_Height;
				currentCol = 0;
			}
			//std::cout << "Done" << std::endl;
			delete[] currentBits;
		}

		for (int i = 0; i < Block_Height; ++i) {
			delete[] block[i];
		}
		delete[] block;
		// Khởi tạo seed dựa trên thời gian hiện tại
		std::srand(std::time(0));

		// Tạo số random từ 0 đến 9
		int random_number = std::rand() % 10;

		// In kết quả
	/*	std::cout << "Free block " << random_number << std::endl;
		std::cout << "Free block hehe " << std::endl;*/

	}

	delete[] encodedBits;

}
