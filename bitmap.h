#ifndef BITMAP_H
#define BITMAP_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

const int FILEHEADERSIZE = 14;
const int INFOHEADERSIZE = 40;
const int HEADERSIZE = 54;

typedef struct{
	unsigned char r,g,b;
} RGB;

class Bitmap {
private:
	unsigned width,height;
	RGB* data;
public:
	Bitmap(int width,int height){
		this->width = width;
		this->height = height;
		data = new RGB[width*height];
		memset(data,0,width*height*sizeof(RGB));
	}
	~Bitmap(){
		delete[] data;
	}
	void set(unsigned x,unsigned y,RGB rgb){
		if(x >= width || y >= height) return;
		data[width * x + y] = rgb;
	}
	void set1(unsigned x,unsigned y){
		RGB rgb;
		rgb.r = rgb.g = rgb.b = (unsigned char)~0u;
		set(x,y,rgb);
	}
	void set0(unsigned x,unsigned y){
		RGB rgb;
		rgb.r = rgb.g = rgb.b = (unsigned char)0u;
		set(x,y,rgb);
	}
	RGB get(unsigned x,unsigned y){
		if(x >= width || y >= height) return data[0];
		return data[width * x + y];	
	}
	void reverse(){
		for(int i = 0; i < (int)(width*height); i++){
			data[i].r = UCHAR_MAX - data[i].r;
			data[i].g = UCHAR_MAX - data[i].g;
			data[i].b = UCHAR_MAX - data[i].b;
		}
	}
	void write(const char* filename){
		int real_width;
		unsigned char* line_data;
		unsigned char header_buf[HEADERSIZE];
		unsigned long file_size;
		unsigned long offset_to_data;
		unsigned long info_header_size;
		unsigned int planes;
		unsigned int color;
		unsigned long compress;
		unsigned data_size;
		long xppm;
		long yppm;

		FILE* fp = fopen(filename,"wb");
		if(fp == NULL){
			fprintf(stderr,"file open error\n");
			return;
		}
		real_width = width*3 + width%4;
		file_size = height * real_width + HEADERSIZE;
		offset_to_data = HEADERSIZE;
		info_header_size = INFOHEADERSIZE;
		planes = 1;
		color = 24;
		compress = 0;
		data_size = height * real_width;
		xppm = 1;
		yppm = 1;

		header_buf[0] = 'B';
		header_buf[1] = 'M';
		memcpy(header_buf+2,&file_size,sizeof(file_size));
		header_buf[6] = 0;
		header_buf[7] = 0;
		header_buf[8] = 0;
		header_buf[9] = 0;
		memcpy(header_buf+10,&offset_to_data,sizeof(offset_to_data));
		memcpy(header_buf+14,&info_header_size,sizeof(info_header_size));
		memcpy(header_buf+18,&width,sizeof(width));
		memcpy(header_buf+22,&height,sizeof(height));
		memcpy(header_buf+26,&planes,sizeof(planes));
		memcpy(header_buf+28,&color,sizeof(color));
		memcpy(header_buf+30,&compress,sizeof(compress));
		memcpy(header_buf+34,&data_size,sizeof(data_size));
		memcpy(header_buf+38,&xppm,sizeof(xppm));
		memcpy(header_buf+42,&yppm,sizeof(yppm));
		memset(header_buf+46,0,8);

		fwrite(header_buf,sizeof(unsigned char),HEADERSIZE,fp);

		line_data = new unsigned char[real_width];
		
		for(int i = 0; i < (int)height; i++){
			for(int j = 0; j < (int)width; j++){
				line_data[j*3]   = data[(height-i-1)*width + j].b;
				line_data[j*3+1] = data[(height-i-1)*width + j].g;
				line_data[j*3+2] = data[(height-i-1)*width + j].r;
			}
			for(int j = width*3;j < real_width; j++){
				line_data[j] = 0;
			}
			fwrite(line_data,sizeof(unsigned char),real_width,fp);
		}
		fclose(fp);
	}
};

#endif
