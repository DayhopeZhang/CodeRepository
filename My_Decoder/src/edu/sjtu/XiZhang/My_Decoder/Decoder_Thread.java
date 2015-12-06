package edu.sjtu.XiZhang.My_Decoder;

import java.util.Arrays;
import android.graphics.Bitmap;
import android.graphics.Bitmap.Config;
import android.graphics.Color;
import android.os.Bundle;
import android.os.Handler;
import android.os.Message;
import android.util.Log;
import com.google.zxing.common.reedsolomon.GenericGF;
import com.google.zxing.common.reedsolomon.ReedSolomonDecoder;
import com.google.zxing.common.reedsolomon.ReedSolomonException;

public class Decoder_Thread extends Thread{
	private char[] rdata = null;
	private static int imgWidth=2448;
	private static int imgHeight=3264;
	final int framSize = imgWidth * imgHeight;
	private Handler uihandler;
	private int[] Data = null;
	
	public Decoder_Thread(int[] data, Handler uihandler){
		this.Data = data;
		this.uihandler = uihandler;
	}
	
	public void run(){
		int[] tmp1 = null,tmp2 = new int[framSize];
		Log.i("TEST", "Starting Decoding.....");
		tmp1 = Binarizer(Data);
	/*	for(int i=0;i<framSize;++i){
			if(tmp1[i] == 0) tmp1[i] = Color.rgb(0, 0, 0);
			else tmp1[i] = Color.rgb(255, 255, 255);
		}*/
		//Log.i("tmp1", String.valueOf(tmp1.length));
		//writeImageToDisk(tmp1,"Binarizer");
		Log.i("TEST", "Binarizer OK!");
		for(int i=0;i<framSize;++i){
			if(tmp1[i] == Color.rgb(0, 0, 0)) tmp1[i] = 0;
			else tmp1[i] = 1;
		}
		dilate(tmp1,tmp2,1);
		Log.i("Dilate", "Dilate OK!");
		erode(tmp2,tmp1,1);
		Log.i("Erode", "Erode OK!");
		/*for(int i=0;i<framSize;++i){
			if(tmp1[i] == 0) tmp1[i] = Color.rgb(0, 0, 0);
			else tmp1[i] = Color.rgb(255, 255, 255);
		}
		writeImageToDisk(tmp1,"Eroded");
		for(int i=0;i<framSize;++i){
			if(tmp1[i] == Color.rgb(0, 0, 0)) tmp1[i] = 0;
			else tmp1[i] = 1;
		}*/
		Point[] p = null;
		p = findCoordinate(tmp1);
		for(int i = 0;i<p.length;++i){
			String str = String.format("X:%d,Y:%d",p[i].x,p[i].y);
			Log.i("Point", str);
		}
		Point[][] p1 = null;
		p1 = findTimingRef(p,tmp1);
		tmp2 = gridOutline(p1,p,tmp1);
		int [] tmp3;
		tmp3 = res_to_decode(tmp2);
		rdata = decode(tmp3);
		
		if(rdata != null){
			Message msg = new Message();
			msg.what = 0x123;
			Bundle mBundle = new Bundle();
			mBundle.putCharArray("rdata", rdata);
			Log.i("returnBytes", Arrays.toString(rdata));
			msg.setData(mBundle);
			uihandler.sendMessage(msg);
		}
	}
	
	/*public static void writeImageToDisk(int[] img, String fileName){
		Bitmap temp = Bitmap.createBitmap(imgWidth, imgHeight, Bitmap.Config.ARGB_8888);
		temp.setPixels(img, 0, imgWidth, 0, 0, imgWidth, imgHeight);
	
		File file = new File(Environment.getExternalStorageDirectory().getPath()+"/DCIM/"+fileName+".JPEG");
		if(!file.exists())
			try {
				file.createNewFile();
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		try{
		    FileOutputStream fops = new FileOutputStream(file);
			temp.compress(Bitmap.CompressFormat.PNG, 100, fops);
			fops.flush();
			fops.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		
		//Toast.makeText(context, "The image stored!", Toast.LENGTH_SHORT).show();
	}*/
	
	
	private int cap(int value,int min,int max){
		if(value<min) return min;
		else if(value>max) return max;
		else return value;
	}
	

	public int[] Binarizer(int[] Y) {
		Bitmap img = Bitmap.createBitmap(Y, imgWidth, imgHeight, Config.ARGB_8888);
		int area = imgWidth * imgHeight;
		int gray[][] = new int[imgWidth][imgHeight];
		int average = 0;// 灰度平均值
		int graysum = 0;
		int graymean = 0;
		int grayfrontmean = 0;
		int graybackmean = 0;
		int pixelGray;
		int front = 0;
		int back = 0;
		int[] pix = new int[imgWidth * imgHeight];
		img.getPixels(pix, 0, imgWidth, 0, 0, imgWidth, imgHeight);
		for (int i = 1; i < imgWidth; i++) { // 不算边界行和列，为避免越界
			for (int j = 1; j < imgHeight; j++) {
				int x = j * imgWidth + i;
				int r = (pix[x] >> 16) & 0xff;
				int g = (pix[x] >> 8) & 0xff;
				int b = pix[x] & 0xff;
				pixelGray = (int) (0.3 * r + 0.59 * g + 0.11 * b);// 计算每个坐标点的灰度
				gray[i][j] = (pixelGray << 16) + (pixelGray << 8) + (pixelGray);
				graysum += pixelGray;
			}
		}
		graymean = (int) (graysum / area);// 整个图的灰度平均值
		average = graymean;
		for (int i = 0; i < imgWidth; i++) // 计算整个图的二值化阈值
		{
			for (int j = 0; j < imgHeight; j++) {
				if (((gray[i][j]) & (0x0000ff)) < graymean) {
					graybackmean += ((gray[i][j]) & (0x0000ff));
					back++;
				} else {
					grayfrontmean += ((gray[i][j]) & (0x0000ff));
					front++;
				}
			}
		}
		int frontvalue = (int) (grayfrontmean / front);// 前景中心
		int backvalue = (int) (graybackmean / back);// 背景中心
		float G[] = new float[frontvalue - backvalue + 1];// 方差数组
		int s = 0;
		
		for (int i1 = backvalue; i1 < frontvalue + 1; i1++)// 以前景中心和背景中心为区间采用大津法算法（OTSU算法）
		{
			back = 0;
			front = 0;
			grayfrontmean = 0;
			graybackmean = 0;
			for (int i = 0; i < imgWidth; i++) {
				for (int j = 0; j < imgHeight; j++) {
					if (((gray[i][j]) & (0x0000ff)) < (i1 + 1)) {
						graybackmean += ((gray[i][j]) & (0x0000ff));
						back++;
					} else {
						grayfrontmean += ((gray[i][j]) & (0x0000ff));
						front++;
					}
				}
			}
			grayfrontmean = (int) (grayfrontmean / front);
			graybackmean = (int) (graybackmean / back);
			G[s] = (((float) back / area) * (graybackmean - average)
					* (graybackmean - average) + ((float) front / area)
					* (grayfrontmean - average) * (grayfrontmean - average));
			s++;
		}
		float max = G[0];
		int index = 0;
		for (int i = 1; i < frontvalue - backvalue + 1; i++) {
			if (max < G[i]) {
				max = G[i];
				index = i;
			}
		}

		for (int i = 0; i < imgWidth; i++) {
			for (int j = 0; j < imgHeight; j++) {
				int in = j * imgWidth + i;
				if (((gray[i][j]) & (0x0000ff)) < (index + backvalue)) {
					pix[in] = Color.rgb(0, 0, 0);
				} else {
					pix[in] = Color.rgb(255, 255, 255);
				}
			}
		}
		return pix;
		//Bitmap temp = Bitmap.createBitmap(imgWidth, imgHeight, Bitmap.Config.ARGB_8888);
		//temp.setPixels(pix, 0, imgWidth, 0, 0, imgWidth, imgHeight);
		//image.setImageBitmap(temp);
    }

	
	public int[] binarizer(int[] Y){
		for(int i = 0; i<framSize;++i) Y[i] = (int)(Color.red(Y[i])*0.2999 + Color.blue(Y[i])*0.114 + Color.green(Y[i])*0.587);
		Log.i("hello", "Hello");
		for(int i=0;i<framSize;++i) Y[i] &= 0x00ff;
		
		int subWidth = imgWidth / 8;
		int subHight = imgHeight /8;
		
		double[] blackPoints = new double[subWidth*subHight];
		int[] tmp = new int[8*8];
		
		for(int i=0;i<subHight;++i){
			for(int j=0;j<subWidth;++j){
				
				for(int k=0;k<8;++k){
					for(int h=0;h<8;++h){
						tmp[k*8+h] = Y[(8*i+k)*imgWidth+(8*j+h)];
					}
				}
				
				int mina = 255;
				
				for(int k=0;k<8;++k){
					if(tmp[k]<mina) mina = tmp[k];
				}
				
				int maxa = 0;
				
				for(int k=0;k<8;++k) {
					if(tmp[k]>maxa) maxa = tmp[k];
				}
				
				int suma = 0;
				if(maxa-mina>24){
					for(int k=0;k<64;++k){
						suma += tmp[k];
					}
				}
				else{
					for(int k=0;k<8;++k){
						suma += tmp[k];
					}
				}
				double average;
				if(maxa-mina <= 24){
					average = mina/2.0;
					
					if(i>0 && j>0){
						double averageNeighborBlackPoint = (blackPoints[(i-1)*subWidth+j] + 2*blackPoints[i*subWidth+j-1] + blackPoints[(i-1)*subWidth+j-1])/4;
						if(mina<averageNeighborBlackPoint) average=averageNeighborBlackPoint;
					}
				}
				else{
					average=suma/64.0;
				}
				blackPoints[i*subWidth+j] = average;
			}
		}
		int[] rdata=new int[imgHeight*imgWidth];
		for(int i=0;i<subHight;++i){
			for(int j=0;j<subWidth;++j){
				int left = cap(j,1,subWidth-4);
				int top = cap(i,1,subHight-4);
				double sum = 0;
				double[] blackRow=new double[subWidth];
				for(int z=-2;z<=2;++z){
					for(int t=0;t<subWidth;++t) blackRow[t] = blackPoints[(top+z+1)*subWidth+t];
					sum = sum + blackRow[left - 2 + 1] + blackRow[left - 1+ 1] + blackRow[left+ 1] + blackRow[left + 1+ 1] + blackRow[left + 2+ 1];
				}
				double average=sum/25;
				
				for(int ii=0;ii<8;++ii){
					for(int jj=0;jj<8;++jj){
						if(Y[(8*i+ii)*imgWidth+8*j+jj] <= average) rdata[(8*i+ii)*imgWidth+8*j+jj] = 1;
					}
				}
			}
		}
		
		for(int i=0;i<imgHeight*imgWidth;++i) rdata[i] = (int)(1-rdata[i]);
		
		return rdata;
	}
	private boolean judge(int m) {if(m>0) return true; else return false;}
	
	private boolean containsBlackPoint(int a,int b, int fixed, boolean h, int[] data){
		//String str = String.format("%d",fixed);
		//Log.i("fixed", str);
		boolean tmp = true;
		boolean flag= false;
		if(h){
			for(int x=a;x<=b;++x){
				tmp = true;
				for(int i=-3;i<=3;++i){
					tmp &= judge(data[x*imgWidth+fixed+i]);
				}
				if(!tmp) {flag=true;break;}
				else flag=false;
			}
		}
		else{
			for(int x=a;x<=b;++x){
				tmp = true;
				for(int i=-3;i<=3;++i){
					tmp &= judge(data[(fixed+i)*imgWidth+x]);
				}
				if(!tmp) {flag = true;break;}
				else flag=false;
			}
		}
		return flag;
	}
	
	private Point getBlackPointOnSegment(int aX,int aY,int bX,int bY,int[] data){
		int x=0,y=0;
		int dis = (int)(Math.sqrt(((aX-bX)*(aX-bX)+(aY-bY)*(aY-bY)))+0.5);
		double xStep = (bX-aX)*1.0 / dis;
		double yStep = (bY-aY)*1.0 / dis;
		for(int i=0;i<=dis;++i){
			x = (int)(aX+i*xStep+0.5);
			y = (int)(aY+i*yStep+0.5);
			if(data[x*imgWidth+y] == 0){
				break;
			}
			else {
				x=0;
				y=0;
			}
		}
		return new Point(x,y);
	}
	
	public Point[] findCoordinate(int[] bitmatrix){
		int left = imgHeight/2-21;
		int right = imgHeight/2+19;
		int up = imgWidth/2-21;
		int down = imgWidth/2+19;
		
		boolean sizeExceed = false,aBlackPointFoundOnBorder = true,atLeastOneBlackPointFoundOnBorder = false,atLeastOneBlackPointFoundOnRight = false;
		boolean atLeastOneBlackPointFoundOnBottom = false,atLeastOneBlackPointFoundOnLeft = false,atLeastOneBlackPointFoundOnTop = false;
		
		while(aBlackPointFoundOnBorder){
			 aBlackPointFoundOnBorder = false;
             boolean rightBorderNotWhite = true;
			 
			 while ((rightBorderNotWhite || !atLeastOneBlackPointFoundOnRight) && right < imgHeight){
				 rightBorderNotWhite = containsBlackPoint(up, down, right,false,bitmatrix);
				 if (rightBorderNotWhite){
					 right=right+1;
                     aBlackPointFoundOnBorder = true;
                     atLeastOneBlackPointFoundOnRight = true;
				 }
				 else if (!atLeastOneBlackPointFoundOnRight){
					 right += 1;
				 }
			 }
			 
			 if (right >= imgHeight){
				 sizeExceed = true;
			     break;
		     }

			 boolean bottomBorderNotWhite = true;
			 while ((bottomBorderNotWhite || !atLeastOneBlackPointFoundOnBottom) && (down < imgWidth)){
				 bottomBorderNotWhite = containsBlackPoint(left, right, down, true,bitmatrix);
				 if (bottomBorderNotWhite){
					 down=down+1;
                     aBlackPointFoundOnBorder = true;
                     atLeastOneBlackPointFoundOnBottom = true;
				 }
				 else if (!atLeastOneBlackPointFoundOnBottom){
					 down = down+1;
				 }
			 }
			 if (down >= imgWidth){
				 sizeExceed=true;
				 break;
			 }
			 
			 boolean leftBorderNotWhite = true;
             while ((leftBorderNotWhite || !atLeastOneBlackPointFoundOnLeft) && left >= 0){
				 leftBorderNotWhite = containsBlackPoint(up, down, left, false,bitmatrix);
				 if (leftBorderNotWhite){
					 left=left-1;
                     aBlackPointFoundOnBorder = true;
                     atLeastOneBlackPointFoundOnLeft = true;
				 }
				 
				 else if(!atLeastOneBlackPointFoundOnLeft){
					 left -= 1;
				 }
			 }
			 if (left < 0){
				 sizeExceed = true;
				 break;
			 }
			 boolean topBorderNotWhite = true;
			 while ((topBorderNotWhite  || !atLeastOneBlackPointFoundOnTop) && up >= 0){
				 topBorderNotWhite = containsBlackPoint(left, right, up, true,bitmatrix);
				 if (topBorderNotWhite){
					 up = up-1;
                     aBlackPointFoundOnBorder = true;
                     atLeastOneBlackPointFoundOnTop = true;
				 }
				 else if (!atLeastOneBlackPointFoundOnTop){
					 up = up-1;
				 }
			 }
			 
			 if(up < 0){
				 sizeExceed=true;
				 break;
			 }
			 if (aBlackPointFoundOnBorder){
				 atLeastOneBlackPointFoundOnBorder = true;
			 }
			 
		}
		
		Point z = new Point(0,0);
		Point t = new Point(0,0);
		Point x = new Point(0,0);
		Point y = new Point(0,0);
		
		String str = String.format("%d",left);
		Log.i("left", str);
		str = String.format("%d",right);
		Log.i("right", str);
		str = String.format("%d",down);
		Log.i("down", str);
		str = String.format("%d",up);
		Log.i("up", str);
		
		if(!sizeExceed && atLeastOneBlackPointFoundOnBorder){
			int maxSize = right - left;
			for(int i = 1;i <= maxSize;++i){
				z = getBlackPointOnSegment(left, down - i, left + i, down,bitmatrix);
				if(z.x != 0 && z.y != 0) break;
			}
			
			for(int i = 1;i <= maxSize;++i){
				t = getBlackPointOnSegment(left, up + i, left + i, up,bitmatrix);
				if(t.x != 0 && t.y != 0) break;
			}
			
			for(int i = 1;i <= maxSize;++i){
				x = getBlackPointOnSegment(right, up + i, right - i, up,bitmatrix);
				if(x.x!=0 && x.y!=0) break;
			}
			
			for(int i = 1;i <= maxSize;++i){
				y = getBlackPointOnSegment(right, down - i, right - i, down,bitmatrix);
				if(y.x!=0 && y.y!=0) break;
			}
		}
		
		double slope = (y.y-t.y)*1.0/(y.x-t.x);
		double revisedis = 2*(Math.sqrt((z.x-t.x)*(z.x-t.x)+(z.y-t.y)*(z.y-t.y))/128/2)*(Math.sqrt((z.x-t.x)*(z.x-t.x)+(z.y-t.y)*(z.y-t.y))/128/2);
		double slope1=(slope-Math.tan(Math.PI/8)) / (Math.tan(Math.PI/8)*slope+1);
		double slope2=(slope+Math.tan(Math.PI/8)) / (-Math.tan(Math.PI/8)*slope+1);
		double A1 = slope1;int B1 = -1;double C1 = z.y-slope1*z.x;
        double A2 = slope2;int B2 = -1;double C2 = z.y-slope2*z.x;
        double A3 = slope1;int B3 = -1;double C3 = x.y-slope1*x.x;
        double A4 = slope2;int B4 = -1;double C4 = x.y-slope2*x.x;
        int m1 = z.x;int n1 = z.y;
        int m2 = z.x;int n2 = z.y;
        int m3 = x.x;int n3 = x.y;
        int m4 = x.x;int n4 = x.y;
        int count1=0,count2=0,count3=0,count4=0;
		
		for(int i=0;i<=(int)(Math.sqrt((z.x-t.x)*(z.x-t.x)+(z.y-t.y)*(z.y-t.y))/130/Math.sqrt(2)*4+0.5);++i){
			if(Math.abs(A1*(m1+1)+B1*n1+C1)/(Math.sqrt(A1*A1+B1*B1))>Math.abs(A1*m1+B1*(n1+1)+C1)/(Math.sqrt(A1*A1+B1*B1))){
				n1 += 1;
			}
			else m1 += 1;
			count1 = count1 + 1 - bitmatrix[m1*imgWidth+n1];  //~(bitmatrix[m1*imgWidth+n1]&1);
			if(Math.abs(A2*(m2-1)+B2*n2+C2)/(Math.sqrt(A2*A2+B2*B2))>Math.abs(A2*m2+B2*(n2-1)+C2)/(Math.sqrt(A2*A2+B2*B2))) n2 -= 1;
			else m2 -= 1;
			count2 = count2 + 1 - bitmatrix[m2*imgWidth+n2]; //~(bitmatrix(m2*imgWidth+n2)&1);
		}
		double zXnew,zYnew;
		if(count1>count2){
			zXnew = Math.sqrt(revisedis/(slope*slope+1))+z.x;
            zYnew = slope*(zXnew-z.x)+z.y;
		}
		else{
			zXnew = -Math.sqrt(revisedis/(slope*slope+1))+z.x;
            zYnew = slope*(zXnew-z.x)+z.y;
		}
		
		for(int i=0;i<=(int)(Math.sqrt((z.x-t.x)*(z.x-t.x)+(z.y-t.y)*(z.y-t.y))/130/Math.sqrt(2)*4+0.5);++i){
			if(Math.abs(A3*(m3-1)+B3*n3+C3)/(Math.sqrt(A3*A3+B3*B3))>Math.abs(A3*m3+B3*(n3-1)+C3)/(Math.sqrt(A3*A3+B3*B3))) n3 -= 1;
			else m3 -= 1;
			count3 = count3 + 1 - bitmatrix[m3*imgWidth+n3]; //~(bitmatrix[m3*imgWidth+n3]&1);
			if(Math.abs(A4*(m4+1)+B4*n4+C4)/(Math.sqrt(A4*A4+B4*B4))>Math.abs(A4*m4+B4*(n4+1)+C4)/(Math.sqrt(A4*A4+B4*B4))) n4 += 1;
			else m4 += 1;
			count4 = count4 + 1 - bitmatrix[m4*imgWidth+n4]; //~(bitmatrix[m4*imgWidth+n4]&1);
		}
		
		double xXnew,xYnew;
		if(count3<count4){
			xXnew = Math.sqrt(revisedis/(slope*slope+1))+x.x;
			xYnew = slope*(xXnew-x.x)+x.y;
		}
		else{
			xXnew = -Math.sqrt(revisedis/(slope*slope+1))+x.x;
            xYnew = slope*(xXnew-x.x)+x.y;
		}
		
		double tXnew = Math.sqrt(revisedis/(slope*slope+1))+t.x;
		double tYnew = slope*(tXnew-t.x)+t.y;
		
		double yXnew = -Math.sqrt(revisedis/(slope*slope+1))+y.x;
        double yYnew = slope*(yXnew-y.x)+y.y;
		
		Point[] rdata = new Point[4];
		for(int i=0; i<4;++i) rdata[i] = new Point();
		
		rdata[0].x = (int)(zXnew+0.5);
		rdata[0].y = (int)(zYnew+0.5);
		rdata[1].x = (int)(tXnew+0.5);
		rdata[1].y = (int)(tYnew+0.5);
		rdata[2].x = (int)(xXnew+0.5);
		rdata[2].y = (int)(xYnew+0.5);
		rdata[3].x = (int)(yXnew+0.5);
		rdata[3].y = (int)(yYnew+0.5);

		return rdata;
	}
	
	public Point[][] findTimingRef(Point[] p,int[] I){
		int black=0,white=1;
		Point pLA = p[1];
		Point pLB = p[0];
		Point pRA = p[2];
		Point pRB = p[3];
		
		Point[] centersL = new Point[67];
		Point[] centersR = new Point[67];
		Point[] centersA = new Point[67];
		Point[] centersB = new Point[67];
		
		for(int i =0; i<67;++i){
			centersL[i] = new Point();
			centersR[i] = new Point();
			centersA[i] = new Point();
			centersB[i] = new Point();
		}
		
		double lenL = Math.sqrt((pLA.x-pLB.x)*(pLA.x-pLB.x)+(pLA.y-pLB.y)*(pLA.y-pLB.y));
		double lenR = Math.sqrt((pRA.x-pRB.x)*(pRA.x-pRB.x)+(pRA.y-pRB.y)*(pRA.y-pRB.y));
		double lenA = Math.sqrt((pLA.x-pRA.x)*(pLA.x-pRA.x)+(pLA.y-pRA.y)*(pLA.y-pRA.y));
		double lenB = Math.sqrt((pLB.x-pRB.x)*(pLB.x-pRB.x)+(pLB.y-pRB.y)*(pLB.y-pRB.y));
		double len = (lenA+lenB+lenL+lenR) / 4;
		
		double sizeEst = len/133;
		int sizeEstInt = (int)(sizeEst+0.5);
		int count=-1;
		
		//Left Pattern
		Point pL = new Point();
		Point p1 = new Point();
		Point p2 = new Point();
		Point p3 = new Point();
		Point p4 = new Point();
		pL.x = pLA.x-pLB.x;
		pL.y = pLA.y-pLB.y;
		int edgetype;
		while(count<66){
			if(count == -1){
				for(int i=1;i<=sizeEstInt;++i){
					p1.x = (int)(pLB.x + pL.x/lenL*i+0.5);
					p1.y = (int)(pLB.y + pL.y/lenL*i+0.5);
					if(I[p1.x*imgWidth+p1.y] == black) break;
				}
			}
			else{
				p1.x = (int)(centersL[count].x + pL.x/lenL*sizeEst+0.5);
				p1.y = (int)(centersL[count].y + pL.y/lenL*sizeEst+0.5);
				Point pBase = new Point();
				pBase.x=p1.x;pBase.y=p1.y;
				for(int i=1;i<=(int)(sizeEstInt*1.5+0.5);++i){
					p1.x = (int)(pBase.x + pL.x/lenL*i+0.5);
					p1.y = (int)(pBase.y + pL.y/lenL*i+0.5);
					if(I[p1.x*imgWidth+p1.y] == black) break;
				}
			}
			
			for(int i=1;i<=(int)(sizeEstInt*1.5+0.5);++i){
				p2.x = (int)(p1.x + pL.x/lenL*i+0.5);
				p2.y = (int)(p1.y + pL.y/lenL*i+0.5);
				if(I[p2.x*imgWidth+p2.y] == white) break;
			}
			p2.x = (int)(p2.x - pL.x/lenL+0.5);
			p2.y = (int)(p2.y - pL.y/lenL+0.5);
			
			Point mi = new Point();
			mi.x = (p1.x+p2.x)/2;
			mi.y = (p1.y+p2.y)/2;
			
			for(int i=1;i<=sizeEstInt;++i){
				p3.x = (int)(mi.x - (-pL.y)/lenL*i+0.5);
				p3.y = (int)(mi.y - (pL.x)/lenL*i+0.5);
				if(I[p3.x*imgWidth+p3.y] == white) break;
			}
			p3.x = (int)(p3.x + (-pL.y)/lenL+0.5);
			p3.y = (int)(p3.y + (pL.x)/lenL+0.5);
			
			edgetype = 0;
			for(int i=1;i<=(int)(sizeEstInt*1.5+0.5);++i){
				p4.x = (int)(p3.x + (-pL.y)/lenL*i+0.5);
				p4.y = (int)(p3.y + (pL.x)/lenL*i+0.5);
				if(I[p4.x*imgWidth+p4.y] == white) {edgetype=1;break;}
			}
			p4.x = (int)(p4.x - (-pL.y)/lenL+0.5);
			p4.y = (int)(p4.y - (pL.x)/lenL+0.5);
			
			if(Math.sqrt((p4.x-p3.x)*(p4.x-p3.x)+(p4.y-p3.y)*(p4.y-p3.y)) > sizeEst){
				p4.x = (int)(p3.x + (-pL.y)/lenL*(sizeEst-1)+0.5);
				p4.y = (int)(p3.y + (pL.x)/lenL*(sizeEst-1)+0.5);
			}
			
			Point c = new Point();
			if(edgetype == 1){
				c.x = (p3.x+p4.x)/2;
				c.y = (p3.y+p4.y)/2;
			}
			else{
				c.x = (int)(p3.x + (mi.x - p3.x)/Math.sqrt((mi.x-p3.x)*(mi.x-p3.x)+(mi.y-p3.y)*(mi.y-p3.y))*(sizeEst/2-1)+0.5);
				c.y = (int)(p3.y + (mi.y - p3.y)/Math.sqrt((mi.x-p3.x)*(mi.x-p3.x)+(mi.y-p3.y)*(mi.y-p3.y))*(sizeEst/2-1)+0.5);
			}
			
			count += 1;
			centersL[count].x = c.x; centersL[count].y = c.y;
		}
		
		// right pattern
		Point pR = new Point();
		pR.x = pRA.x-pRB.x;
		pR.y = pRA.y-pRB.y;
		count = -1;
		while(count<66){
			if(count == -1){
				for(int i=1;i<=sizeEstInt;++i){
					p1.x = (int)(pRB.x + pR.x/lenR*i+0.5);
					p1.y = (int)(pRB.y + pR.y/lenR*i+0.5);
					if(I[p1.x*imgWidth+p1.y] == black) break;
				}
			}
			else{
				p1.x = (int)(centersR[count].x + pR.x/lenR*sizeEst+0.5);
				p1.y = (int)(centersR[count].y + pR.y/lenR*sizeEst+0.5);
				Point pBase = new Point();
				pBase.x=p1.x;pBase.y=p1.y;
				for(int i=1;i<=(int)(sizeEstInt*1.5+0.5);++i){
					p1.x = (int)(pBase.x + pR.x/lenR*i+0.5);
					p1.y = (int)(pBase.y + pR.y/lenR*i+0.5);
					if(I[p1.x*imgWidth+p1.y] == black) break;
				}
			}
			
			for(int i=1;i<=(int)(sizeEstInt*1.5+0.5);++i){
				p2.x = (int)(p1.x + pR.x/lenR*i+0.5);
				p2.y = (int)(p1.y + pR.y/lenR*i+0.5);
				if(I[p2.x*imgWidth+p2.y] == white) break;
			}
			p2.x = (int)(p2.x - pR.x/lenR+0.5);
			p2.y = (int)(p2.y - pR.y/lenR+0.5);
			
			Point mi = new Point();
			mi.x = (p1.x+p2.x)/2;
			mi.y = (p1.y+p2.y)/2;
			
			for(int i=1;i<=sizeEstInt;++i){
				p3.x = (int)(mi.x + (-pR.y)/lenR*i+0.5);
				p3.y = (int)(mi.y + (pR.x)/lenR*i+0.5);
				if(I[p3.x*imgWidth+p3.y] == white) break;
			}
			p3.x = (int)(p3.x - (-pR.y)/lenR+0.5);
			p3.y = (int)(p3.y - (pR.x)/lenR+0.5);
			
			edgetype = 0;
			for(int i=1;i<=(int)(sizeEstInt*1.5+0.5);++i){
				p4.x = (int)(p3.x - (-pR.y)/lenR*i+0.5);
				p4.y = (int)(p3.y - (pR.x)/lenR*i+0.5);
				if(I[p4.x*imgWidth+p4.y] == white) {edgetype=1;break;}
			}
			p4.x = (int)(p4.x + (-pR.y)/lenR+0.5);
			p4.y = (int)(p4.y + (pR.x)/lenR+0.5);
			
			if(Math.sqrt((p4.x-p3.x)*(p4.x-p3.x)+(p4.y-p3.y)*(p4.y-p3.y)) > sizeEst){
				p4.x = (int)(p3.x - (-pR.y)/lenR*(sizeEst-1)+0.5);
				p4.y = (int)(p3.y - (pR.x)/lenR*(sizeEst-1)+0.5);
			}
			
			Point c = new Point();
			if(edgetype == 1){
				c.x = (p3.x+p4.x)/2;
				c.y = (p3.y+p4.y)/2;
			}
			else{
				c.x = (int)(p3.x + (mi.x - p3.x)/Math.sqrt((mi.x-p3.x)*(mi.x-p3.x)+(mi.y-p3.y)*(mi.y-p3.y))*(sizeEst/2-1)+0.5);
				c.y = (int)(p3.y + (mi.y - p3.y)/Math.sqrt((mi.x-p3.x)*(mi.x-p3.x)+(mi.y-p3.y)*(mi.y-p3.y))*(sizeEst/2-1)+0.5);
			}
			
			count += 1;
			centersR[count].x = c.x; centersR[count].y = c.y;
			centersR[0].x = pRB.x; centersR[0].y = pRB.y;
		}
		
		//top pattern
		Point pA = new Point();
		pA.x = pRA.x-pLA.x;
		pA.y = pRA.y-pLA.y;
		count = -1;
		while(count<66){
			if(count == -1){
				for(int i=1;i<=sizeEstInt;++i){
					p1.x = (int)(pLA.x + pA.x/lenA*i+0.5);
					p1.y = (int)(pLA.y + pA.y/lenA*i+0.5);
					if(I[p1.x*imgWidth+p1.y] == black) break;
				}
			}
			else{
				p1.x = (int)(centersA[count].x + pA.x/lenA*sizeEst+0.5);
				p1.y = (int)(centersA[count].y + pA.y/lenA*sizeEst+0.5);
				Point pBase = new Point();
				pBase.x=p1.x;pBase.y=p1.y;
				for(int i=1;i<=(int)(sizeEstInt*1.5+0.5);++i){
					p1.x = (int)(pBase.x + pA.x/lenA*i+0.5);
					p1.y = (int)(pBase.y + pA.y/lenA*i+0.5);
					if(I[p1.x*imgWidth+p1.y] == black) break;
				}
			}
			
			for(int i=1;i<=(int)(sizeEstInt*1.5+0.5);++i){
				p2.x = (int)(p1.x + pA.x/lenA*i+0.5);
				p2.y = (int)(p1.y + pA.y/lenA*i+0.5);
				if(I[p2.x*imgWidth+p2.y] == white) break;
			}
			p2.x = (int)(p2.x - pA.x/lenA+0.5);
			p2.y = (int)(p2.y - pA.y/lenA+0.5);
			
			Point mi = new Point();
			mi.x = (p1.x+p2.x)/2;
			mi.y = (p1.y+p2.y)/2;
			
			for(int i=1;i<=sizeEstInt;++i){
				p3.x = (int)(mi.x - (-pA.y)/lenA*i+0.5);
				p3.y = (int)(mi.y - (pA.x)/lenA*i+0.5);
				if(I[p3.x*imgWidth+p3.y] == white) break;
			}
			p3.x = (int)(p3.x + (-pA.y)/lenA+0.5);
			p3.y = (int)(p3.y + (pA.x)/lenA+0.5);
			
			edgetype = 0;
			for(int i=1;i<=(int)(sizeEstInt*1.5+0.5);++i){
				p4.x = (int)(p3.x + (-pA.y)/lenA*i+0.5);
				p4.y = (int)(p3.y + (pA.x)/lenA*i+0.5);
				if(I[p4.x*imgWidth+p4.y] == white) {edgetype=1;break;}
			}
			p4.x = (int)(p4.x - (-pA.y)/lenA+0.5);
			p4.y = (int)(p4.y - (pA.x)/lenA+0.5);
			
			if(Math.sqrt((p4.x-p3.x)*(p4.x-p3.x)+(p4.y-p3.y)*(p4.y-p3.y)) > sizeEst){
				p4.x = (int)(p3.x + (-pA.y)/lenA*(sizeEst-1)+0.5);
				p4.y = (int)(p3.y + (pA.x)/lenA*(sizeEst-1)+0.5);
			}
			
			Point c = new Point();
			if(edgetype == 1){
				c.x = (p3.x+p4.x)/2;
				c.y = (p3.y+p4.y)/2;
			}
			else{
				c.x = (int)(p3.x + (mi.x - p3.x)/Math.sqrt((mi.x-p3.x)*(mi.x-p3.x)+(mi.y-p3.y)*(mi.y-p3.y))*(sizeEst/2-1)+0.5);
				c.y = (int)(p3.y + (mi.y - p3.y)/Math.sqrt((mi.x-p3.x)*(mi.x-p3.x)+(mi.y-p3.y)*(mi.y-p3.y))*(sizeEst/2-1)+0.5);
			}
			
			count += 1;
			centersA[count].x = c.x; centersA[count].y = c.y;
			centersA[0].x = pLA.x; centersA[0].y = pLA.y;
		}
		
		// bottom pattern
		Point pB = new Point();
		pB.x = pRB.x-pLB.x;
		pB.y = pRB.y-pLB.y;
		count = -1;
		while(count<66){
			if(count == -1){
				for(int i=1;i<=sizeEstInt;++i){
					p1.x = (int)(pLB.x + pB.x/lenB*i+0.5);
					p1.y = (int)(pLB.y + pB.y/lenB*i+0.5);
					if(I[p1.x*imgWidth+p1.y] == black) break;
				}
			}
			else{
				p1.x = (int)(centersB[count].x + pB.x/lenB*sizeEst+0.5);
				p1.y = (int)(centersB[count].y + pB.y/lenB*sizeEst+0.5);
				Point pBase = new Point();
				pBase.x=p1.x;pBase.y=p1.y;
				for(int i=1;i<=(int)(sizeEstInt*1.5+0.5);++i){
					p1.x = (int)(pBase.x + pB.x/lenB*i+0.5);
					p1.y = (int)(pBase.y + pB.y/lenB*i+0.5);
					if(I[p1.x*imgWidth+p1.y] == black) break;
				}
			}
			
			for(int i=1;i<=(int)(sizeEstInt*1.5+0.5);++i){
				p2.x = (int)(p1.x + pB.x/lenB*i+0.5);
				p2.y = (int)(p1.y + pB.y/lenB*i+0.5);
				if(I[p2.x*imgWidth+p2.y] == white) break;
			}
			p2.x = (int)(p2.x - pB.x/lenB+0.5);
			p2.y = (int)(p2.y - pB.y/lenB+0.5);
			
			Point mi = new Point();
			mi.x = (p1.x+p2.x)/2;
			mi.y = (p1.y+p2.y)/2;
			
			for(int i=1;i<=sizeEstInt;++i){
				p3.x = (int)(mi.x + (-pB.y)/lenB*i+0.5);
				p3.y = (int)(mi.y + (pB.x)/lenB*i+0.5);
				if(I[p3.x*imgWidth+p3.y] == white) break;
			}
			p3.x = (int)(p3.x - (-pB.y)/lenB+0.5);
			p3.y = (int)(p3.y - (pB.x)/lenB+0.5);
			
			edgetype = 0;
			for(int i=1;i<=(int)(sizeEstInt*1.5+0.5);++i){
				p4.x = (int)(p3.x - (-pB.y)/lenB*i+0.5);
				p4.y = (int)(p3.y - (pB.x)/lenB*i+0.5);
				if(I[p4.x*imgWidth+p4.y] == white) {edgetype=1;break;}
			}
			p4.x = (int)(p4.x + (-pB.y)/lenB+0.5);
			p4.y = (int)(p4.y + (pB.x)/lenB+0.5);
			
			if(Math.sqrt((p4.x-p3.x)*(p4.x-p3.x)+(p4.y-p3.y)*(p4.y-p3.y)) > sizeEst){
				p4.x = (int)(p3.x - (-pB.y)/lenB*(sizeEst-1)+0.5);
				p4.y = (int)(p3.y - (pB.x)/lenB*(sizeEst-1)+0.5);
			}
			
			Point c = new Point();
			if(edgetype == 1){
				c.x = (p3.x+p4.x)/2;
				c.y = (p3.y+p4.y)/2;
			}
			else{
				c.x = (int)(p3.x + (mi.x - p3.x)/Math.sqrt((mi.x-p3.x)*(mi.x-p3.x)+(mi.y-p3.y)*(mi.y-p3.y))*(sizeEst/2-1)+0.5);
				c.y = (int)(p3.y + (mi.y - p3.y)/Math.sqrt((mi.x-p3.x)*(mi.x-p3.x)+(mi.y-p3.y)*(mi.y-p3.y))*(sizeEst/2-1)+0.5);
			}
			
			count += 1;
			centersB[count].x = c.x; centersB[count].y = c.y;
		}
		
		Point[][] rdata = new Point[4][67];
		for(int i=0;i<67;++i){
			rdata[0][i] = centersL[i];
			rdata[1][i] = centersR[i];
			rdata[2][i] = centersA[i];
			rdata[3][i] = centersB[i];
		}
		return rdata;
	}
	
	public int[] gridOutline(Point[][] data,Point[] p, int[] I){
		Point pLA = p[1];
		Point pLB = p[0];
		Point pRA = p[2];
		Point pRB = p[3];
		Point[] centersL = new Point[67];
		Point[] centersR = new Point[67];
		Point[] centersA = new Point[67];
		Point[] centersB = new Point[67];
		int[] M = new int[132*132];
		
		for(int i=0;i<67;++i){
			centersL[i] = data[0][i];
			centersR[i] = data[1][i];
			centersA[i] = data[2][i];
			centersB[i] = data[3][i];
		}
		
		Point[] gridLef = new Point[132];
		Point[] gridRig = new Point[132];
		Point[] gridTop = new Point[132];
		Point[] gridBot = new Point[132];
		
		for(int i = 0;i < 132;++i){
			gridLef[i] = new Point();
			gridRig[i] = new Point();
			gridTop[i] = new Point();
			gridBot[i] = new Point();
		}
		
		for(int i=0;i<66;++i){
			gridLef[2*i] = centersL[i];
			gridLef[2*i+1].x = (int)((centersL[i].x+centersL[i+1].x)/2.0+0.5);
			gridLef[2*i+1].y = (int)((centersL[i].y+centersL[i+1].y)/2.0+0.5);
			gridRig[2*i].x = (int)((centersR[i].x+centersR[i+1].x)/2.0+0.5);
			gridRig[2*i].y = (int)((centersR[i].y+centersR[i+1].y)/2.0+0.5);
			gridRig[2*i+1] = centersR[i+1];
			gridTop[2*i].x = (int)((centersA[i].x+centersA[i+1].x)/2.0+0.5);
			gridTop[2*i].y = (int)((centersA[i].y+centersA[i+1].y)/2.0+0.5);
			gridTop[2*i+1] = centersA[i+1];
			gridBot[2*i] = centersB[i];
			gridBot[2*i+1].x = (int)((centersB[i].x+centersB[i+1].x)/2.0+0.5);
			gridBot[2*i+1].y = (int)((centersB[i].y+centersB[i+1].y)/2.0+0.5);
		}
		
		for(int i=0;i<132;++i){
			for(int j=0;j<132;++j){
				double left = (gridLef[131-j].y - pLA.y)*1.0 / (pLB.y-pLA.y);
				double right = (gridRig[131-j].y - pRA.y)*1.0 / (pRB.y-pRA.y);
				double rx1 = (gridTop[i].x-pLA.x)*1.0 / (pRA.x-pLA.x);
				double ry1 = left * (1-rx1) + right * rx1;
				double y = gridTop[i].y * (1-ry1) + gridBot[i].y * ry1;
				
				double top = (gridTop[i].x - pRA.x)*1.0 / (pLA.x-pRA.x);
                double bottom = (gridBot[i].x-pRB.x)*1.0 / (pLB.x-pRB.x);
                double rx2 = (gridRig[131-j].y-pRA.y)*1.0 / (pRB.y-pRA.y);
                double ry2 = top * (1-rx2) + bottom * rx2;
                double x = gridRig[131-j].x * (1-ry2) + gridLef[131-j].x * ry2;
				
				M[i*132+j] = I[((int)(x+0.5))*imgWidth+(int)(y+0.5)];
			}
		}
		
		return M;
	}
	
	public int[] res_to_decode(int[] M){
		/*int[] M = new int[132*132];
		for(int i=0;i<132;++i){
			for(int j=0;j<132;++j){
				M[i*132+j] = T[j*132+i];
			}
		}*/
		int[] A = new int[256*256];
		
		for(int i = 0;i<128;++i){
			for(int j=0;j<128;++j){
				A[2*i*256+2*j] = M[(i+3)*132+j+2];
				A[2*i*256+(2*j+1)] = M[(i+3)*132+j+2];
				A[(2*i+1)*256+2*j] = M[(i+3)*132+j+2];
				A[(2*i+1)*256+(2*j+1)] = M[(i+3)*132+j+2];
			}
		}
		
		return A;
	}
	
	public int[] rot90(int[] a){
		int[] rdata = new int[6*6];
		for(int i=0;i<6;++i){
			for(int j=0;j<6;++j){
				rdata[(5-j)*6+i] = a[i*6+j];
			}
		}
		return rdata;
	}
	
	public int[] flipud(int[] a){
		int[] rdata = new int[6*6];
		for(int i = 0;i<6;++i){
			for(int j = 0;j<6;++j){
				rdata[i*6+j] = a[(5-i)*6+j];
			}
		}
		return rdata;
	}
	
	public char[] decode(int[] A){
		int[] II = new int[256*256];
		for(int i=0;i<256*256;++i) II[i] = A[i];
		int[] row = new int[1000];
		int[] col = new int[1000];
		int point = 0, d = 0, c=0;
		
		int[] x1={1,1,d,d,1,1,1,1,d,d,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
		int[] x2 = rot90(x1);
		int[] x3 = rot90(x2);
		int[] x4 = rot90(x3);
		
		int[] x5={d,d,1,1,d,d,d,d,1,1,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d};
		int[] x6=rot90(x5);
		int[] x7=rot90(x6);
		int[] x8=rot90(x7);
		
		int[] x9={1,1,d,d,1,1,1,1,d,d,1,1,1,1,d,d,1,1,1,1,d,d,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
		int[] x10 = rot90(x9);
		int[] x11 = rot90(x10);
		int[] x12 = rot90(x11);
		
		int[] x13={d,d,1,1,d,d,d,d,1,1,d,d,d,d,1,1,d,d,d,d,1,1,d,d,d,d,d,d,d,d,d,d,d,d,d,d};
		int[] x14=rot90(x13);
		int[] x15=rot90(x14);
		int[] x16=rot90(x15);
		
		int[] y1={1,1,d,d,d,d,1,1,d,d,d,d,1,1,d,d,d,d,1,1,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d};
		int[] y2 = rot90(y1);
        int[] y3 = rot90(y2);
        int[] y4 = rot90(y3);
        int[] y5 = flipud(y1);
        int[] y6 = rot90(y5);
        int[] y7 = rot90(y6);
        int[] y8 = rot90(y7);
		
		int[] y9={1,1,d,d,d,d,1,1,d,d,d,d,1,1,1,1,d,d,1,1,1,1,d,d,d,d,d,d,d,d,d,d,d,d,d,d};
		int[] y10 = rot90(y9);
        int[] y11 = rot90(y10);
        int[] y12 = rot90(y11);
        int[] y13 = flipud(y9);
        int[] y14 = rot90(y13);
        int[] y15 = rot90(y14);
        int[] y16 = rot90(y15);
		
		int[] z1={d,d,1,1,1,1,d,d,1,1,1,1,d,d,1,1,1,1,d,d,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
		int[] z2 = rot90(z1);
        int[] z3 = rot90(z2);
        int[] z4 = rot90(z3);
        int[] z5 = flipud(z1);
        int[] z6 = rot90(z5);
        int[] z7 = rot90(z6);
        int[] z8 = rot90(z7);
		
		int[] z9={d,d,1,1,1,1,d,d,1,1,1,1,d,d,d,d,1,1,d,d,d,d,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
		int[] z10 = rot90(z9);
        int[] z11 = rot90(z10);
        int[] z12 = rot90(z11);
        int[] z13 = flipud(z9);
        int[] z14 = rot90(z13);
        int[] z15 = rot90(z14);
        int[] z16 = rot90(z15);
		
		int[] a1={d,d,1,1,1,1,d,d,1,1,1,1,1,1,d,d,1,1,1,1,d,d,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
		int[] a2 = rot90(a1);
        int[] a3 = rot90(a2);
        int[] a4 = rot90(a3);
		
		int[] a5={d,d,1,1,1,1,d,d,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
		int[] a6 = rot90(a5);
        int[] a7 = rot90(a6);
        int[] a8 = rot90(a7);
        
        int[] b1={d,d,1,1,1,1,d,d,1,1,1,1,d,d,1,1,1,1,d,d,1,1,1,1,d,d,1,1,1,1,d,d,1,1,1,1};
        int[] b2 = rot90(b1);
        int[] b3 = rot90(b2);
        int[] b4 = rot90(b3);
        
        int[] b5={d,d,1,1,1,1,d,d,1,1,1,1,d,d,d,d,1,1,d,d,d,d,1,1,d,d,1,1,1,1,d,d,1,1,1,1};
        int[] b6 = rot90(b5);
        int[] b7 = rot90(b6);
        int[] b8 = rot90(b7);
        
        int[] c1={d,d,1,1,d,d,d,d,1,1,d,d,d,d,1,1,1,1,d,d,1,1,1,1,1,1,d,d,1,1,1,1,d,d,1,1};
        int[] c2 = rot90(c1);
        int[] c3 = rot90(c2);
        int[] c4 = rot90(c3);
        int[] c5 = flipud(c1);
        int[] c6 = rot90(c5);
        int[] c7 = rot90(c6);
        int[] c8 = rot90(c7);

        int[] c9={d,d,1,1,d,d,d,d,1,1,d,d,d,d,d,d,1,1,d,d,d,d,1,1,1,1,d,d,1,1,1,1,d,d,1,1};
        int[] c10 = rot90(c9);
        int[] c11 = rot90(c10);
        int[] c12 = rot90(c11);
        int[] c13 = flipud(c9);
        int[] c14 = rot90(c13);
        int[] c15 = rot90(c14);
        int[] c16 = rot90(c15);

        int[] d1={d,d,1,1,d,d,d,d,1,1,d,d,1,1,1,1,d,d,1,1,1,1,d,d,1,1,d,d,d,d,1,1,d,d,d,d};
        int[] d2 = rot90(d1);
        int[] d3 = rot90(d2);
        int[] d4 = rot90(d3);
        int[] d5 = flipud(d1);
        int[] d6 = rot90(d5);
        int[] d7 = rot90(d6);
        int[] d8 = rot90(d7);
       
        int[] d9={d,d,1,1,d,d,d,d,1,1,d,d,1,1,d,d,d,d,1,1,d,d,d,d,1,1,d,d,d,d,1,1,d,d,d,d};
        int[] d10 = rot90(d9);
        int[] d11 = rot90(d10);
        int[] d12 = rot90(d11);
        int[] d13 = flipud(d9);
        int[] d14 = rot90(d13);
        int[] d15 = rot90(d14);
        int[] d16 = rot90(d15);

        for(int l=0;l<256;++l){
			for(int k=0;k<256;++k){
				if(A[l*256+k] == 0) II[l*256+k] = d;
				else A[l*256+k] = 1;
			}
		}
        
        int num=0;
        for(int i=0;i<249;i=i+4){
			for(int j=0;j<249;j=j+4){
		    boolean flag1=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != x1[k*6+h]) {flag1=false;break;}
				}
			}
			boolean flag2=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != x2[k*6+h]) {flag2=false;break;}
				}
			}
			boolean flag3=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != x3[k*6+h]) {flag3=false;break;}
				}
			}
			boolean flag4=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != x4[k*6+h]) {flag4=false;break;}
				}
			}
			boolean flag5=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != x5[k*6+h]) {flag5=false;break;}
				}
			}
			boolean flag6=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != x6[k*6+h]) {flag6=false;break;}
				}
			}
			boolean flag7=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != x7[k*6+h]) {flag7=false;break;}
				}
			}
			boolean flag8=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != x8[k*6+h]) {flag8=false;break;}
				}
			}
			boolean flag9=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != x9[k*6+h]) {flag9=false;break;}
				}
			}
			boolean flag10=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != x10[k*6+h]) {flag10=false;break;}
				}
			}
			boolean flag11=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != x11[k*6+h]) {flag11=false;break;}
				}
			}
			boolean flag12=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != x12[k*6+h]) {flag12=false;break;}
				}
			}
			boolean flag13=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != x13[k*6+h]) {flag13=false;break;}
				}
			}
			boolean flag14=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != x14[k*6+h]) {flag14=false;break;}
				}
			}
			boolean flag15=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != x15[k*6+h]) {flag15=false;break;}
				}
			}
			boolean flag16=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != x16[k*6+h]) {flag16=false;break;}
				}
			}
			boolean flag17=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != y1[k*6+h]) {flag17=false;break;}
				}
			}
			boolean flag18=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != y2[k*6+h]) {flag18=false;break;}
				}
			}
			boolean flag19=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != y3[k*6+h]) {flag19=false;break;}
				}
			}
			boolean flag20=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != y4[k*6+h]) {flag20=false;break;}
				}
			}
			boolean flag21=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != y5[k*6+h]) {flag21=false;break;}
				}
			}
			boolean flag22=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != y6[k*6+h]) {flag22=false;break;}
				}
			}
			boolean flag23=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != y7[k*6+h]) {flag23=false;break;}
				}
			}
			boolean flag24=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != y8[k*6+h]) {flag24=false;break;}
				}
			}
			boolean flag25=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != y9[k*6+h]) {flag25=false;break;}
				}
			}
			boolean flag26=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != y10[k*6+h]) {flag26=false;break;}
				}
			}
			boolean flag27=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != y11[k*6+h]) {flag27=false;break;}
				}
			}
			boolean flag28=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != y12[k*6+h]) {flag28=false;break;}
				}
			}
			boolean flag29=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != y13[k*6+h]) {flag29=false;break;}
				}
			}
			boolean flag30=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != y14[k*6+h]) {flag30=false;break;}
				}
			}
			boolean flag31=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != y15[k*6+h]) {flag31=false;break;}
				}
			}
			boolean flag32=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != y16[k*6+h]) {flag32=false;break;}
				}
			}
			boolean flag33=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != z1[k*6+h]) {flag33=false;break;}
				}
			}
			boolean flag34=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != z2[k*6+h]) {flag34=false;break;}
				}
			}
			boolean flag35=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != z3[k*6+h]) {flag35=false;break;}
				}
			}
			boolean flag36=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != z4[k*6+h]) {flag36=false;break;}
				}
			}
			boolean flag37=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != z5[k*6+h]) {flag37=false;break;}
				}
			}
			boolean flag38=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != z6[k*6+h]) {flag38=false;break;}
				}
			}
			boolean flag39=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != z7[k*6+h]) {flag39=false;break;}
				}
			}
			boolean flag40=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != z8[k*6+h]) {flag40=false;break;}
				}
			}
			boolean flag41=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != z9[k*6+h]) {flag41=false;break;}
				}
			}
			boolean flag42=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != z10[k*6+h]) {flag42=false;break;}
				}
			}
			boolean flag43=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != z11[k*6+h]) {flag43=false;break;}
				}
			}
			boolean flag44=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != z12[k*6+h]) {flag44=false;break;}
				}
			}
			boolean flag45=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != z13[k*6+h]) {flag45=false;break;}
				}
			}
			boolean flag46=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != z14[k*6+h]) {flag46=false;break;}
				}
			}
			boolean flag47=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != z15[k*6+h]) {flag47=false;break;}
				}
			}
			boolean flag48=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != z16[k*6+h]) {flag48=false;break;}
				}
			}
			boolean flag49=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != a1[k*6+h]) {flag49=false;break;}
				}
			}
			boolean flag50=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != a2[k*6+h]) {flag50=false;break;}
				}
			}
			boolean flag51=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != a3[k*6+h]) {flag51=false;break;}
				}
			}
			boolean flag52=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != a4[k*6+h]) {flag52=false;break;}
				}
			}
			boolean flag53=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != a5[k*6+h]) {flag53=false;break;}
				}
			}
			boolean flag54=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != a6[k*6+h]) {flag54=false;break;}
				}
			}
			boolean flag55=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != a7[k*6+h]) {flag55=false;break;}
				}
			}
			boolean flag56=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != a8[k*6+h]) {flag56=false;break;}
				}
			}
			boolean flag57=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != b1[k*6+h]) {flag57=false;break;}
				}
			}
			boolean flag58=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != b2[k*6+h]) {flag58=false;break;}
				}
			}
			boolean flag59=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != b3[k*6+h]) {flag59=false;break;}
				}
			}
			boolean flag60=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != b4[k*6+h]) {flag60=false;break;}
				}
			}
			boolean flag61=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != b5[k*6+h]) {flag61=false;break;}
				}
			}
			boolean flag62=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != b6[k*6+h]) {flag62=false;break;}
				}
			}
			boolean flag63=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != b7[k*6+h]) {flag63=false;break;}
				}
			}
			boolean flag64=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != b8[k*6+h]) {flag64=false;break;}
				}
			}
			boolean flag65=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != c1[k*6+h]) {flag65=false;break;}
				}
			}
			boolean flag66=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != c2[k*6+h]) {flag66=false;break;}
				}
			}
			boolean flag67=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != c3[k*6+h]) {flag67=false;break;}
				}
			}
			boolean flag68=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != c4[k*6+h]) {flag68=false;break;}
				}
			}
			boolean flag69=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != c5[k*6+h]) {flag69=false;break;}
				}
			}
			boolean flag70=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != c6[k*6+h]) {flag70=false;break;}
				}
			}
			boolean flag71=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != c7[k*6+h]) {flag71=false;break;}
				}
			}
			boolean flag72=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != c8[k*6+h]) {flag72=false;break;}
				}
			}
			boolean flag73=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != c9[k*6+h]) {flag73=false;break;}
				}
			}
			boolean flag74=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != c10[k*6+h]) {flag74=false;break;}
				}
			}
			boolean flag75=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != c11[k*6+h]) {flag75=false;break;}
				}
			}
			boolean flag76=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != c12[k*6+h]) {flag76=false;break;}
				}
			}
			boolean flag77=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != c13[k*6+h]) {flag77=false;break;}
				}
			}
			boolean flag78=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != c14[k*6+h]) {flag78=false;break;}
				}
			}
			boolean flag79=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != c15[k*6+h]) {flag79=false;break;}
				}
			}
			boolean flag80=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != c16[k*6+h]) {flag80=false;break;}
				}
			}
			boolean flag81=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != d1[k*6+h]) {flag81=false;break;}
				}
			}
			boolean flag82=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != d2[k*6+h]) {flag82=false;break;}
				}
			}
			boolean flag83=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != d3[k*6+h]) {flag83=false;break;}
				}
			}
			boolean flag84=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != d4[k*6+h]) {flag84=false;break;}
				}
			}
			boolean flag85=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != d5[k*6+h]) {flag85=false;break;}
				}
			}
			boolean flag86=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != d6[k*6+h]) {flag86=false;break;}
				}
			}
			boolean flag87=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != d7[k*6+h]) {flag87=false;break;}
				}
			}
			boolean flag88=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != d8[k*6+h]) {flag88=false;break;}
				}
			}
			boolean flag89=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != d9[k*6+h]) {flag89=false;break;}
				}
			}
			boolean flag90=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != d10[k*6+h]) {flag90=false;break;}
				}
			}
			boolean flag91=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != d11[k*6+h]) {flag91=false;break;}
				}
			}
			boolean flag92=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != d12[k*6+h]) {flag92=false;break;}
				}
			}
			boolean flag93=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != d13[k*6+h]) {flag93=false;break;}
				}
			}
			boolean flag94=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != d14[k*6+h]) {flag94=false;break;}
				}
			}
			boolean flag95=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != d15[k*6+h]) {flag95=false;break;}
				}
			}
			boolean flag96=true;
			for(int k=0;k<6;++k){
				for(int h=0;h<6;++h){
					if(II[(i+k)*256+j+h] != d16[k*6+h]) {flag96=false;break;}
				}
			}
			if(II[(i+2)*256+j]==d && II[(i+3)*256+j]==d && II[(i+2)*256+j+1] == d && II[(i+3)*256+j+1] == d && II[i*256+j+2] == d && II[i*256+j+3] == d && II[(i+1)*256+j+2] == d && II[(i+1)*256+j+3] == d && II[(i+2)*256+j+4] == 1 && II[(i+2)*256+j+5] == 1 && II[(i+3)*256+j+4] == 1 && II[(i+3)*256+j+5] == 1 && II[(i+4)*256+j+2] == 1 && II[(i+4)*256+j+3] == 1 && II[(i+5)*256+j+2] == 1 && II[(i+5)*256+j+3] == 1 ){
				num = num + 1;  
                row[point] = i + 2;
                col[point] = j + 2;
                point = point + 1;
			}
			else if(II[(i+2)*256+j] == d && II[(i+3)*256+j] == d && II[(i+2)*256+j+1] == d && II[(i+3)*256+j+1] == d && II[(i+4)*256+j+2] == d && II[(i+4)*256+j+3] == d && II[(i+5)*256+j+2] == d && II[(i+5)*256+j+3] == d && II[(i+2)*256+j+4] == 1 && II[(i+2)*256+j+5] == 1 && II[(i+3)*256+j+4] == 1 && II[(i+3)*256+j+5] == 1 && II[i*256+j+2] == 1 && II[i*256+j+3] == 1 && II[(i+1)*256+j+2] == 1 && II[(i+1)*256+j+3] == 1 ){
				num = num + 1;
                row[point] = i + 2;
                col[point] = j + 2;
                point = point + 1;
			}
			else if(II[(i+2)*256+j+4] == d && II[(i+3)*256+j+4] == d && II[(i+2)*256+j+5] == d && II[(i+3)*256+j+5] == d && II[(i+4)*256+j+2] == d && II[(i+4)*256+j+3] == d && II[(i+5)*256+j+2] == d && II[(i+5)*256+j+3] == d && II[(i+2)*256+j] == 1 && II[(i+2)*256+j+1] == 1 && II[(i+3)*256+j] == 1 && II[(i+3)*256+j+1] == 1 && II[i*256+j+2] == 1 && II[i*256+j+3] == 1 && II[(i+1)*256+j+2] == 1 && II[(i+1)*256+j+3] == 1){
				num = num + 1;
                row[point] = i + 2;
                col[point] = j + 2;
                point = point + 1;
			}
			else if(II[(i+2)*256+j+4] == d && II[(i+3)*256+j+4] == d && II[(i+2)*256+j+5] == d && II[(i+3)*256+j+5] == d && II[i*256+j+2] == d && II[i*256+j+3] == d && II[(i+1)*256+j+2] == d && II[(i+1)*256+j+3] == d && II[(i+2)*256+j] == 1 && II[(i+2)*256+j+1] == 1 && II[(i+3)*256+j] == 1 && II[(i+3)*256+j+1] == 1 && II[(i+4)*256+j+2] == 1 && II[(i+4)*256+j+3] == 1 && II[(i+5)*256+j+2] == 1 && II[(i+5)*256+j+3] == 1){
				num = num + 1;
                row[point] = i + 2;
                col[point] = j + 2;
                point = point + 1;
			}
			else if(flag1 || flag2 || flag3 || flag4 || flag5 || flag6 || flag7 || flag8 || flag9 || flag10 || flag11 || flag12 || flag13 || flag14 || flag15 || flag16 || flag17 || flag18 || flag19 || flag20 || flag21 || flag22 || flag23 || flag24 || flag25 || flag26 || flag27 || flag28 || flag29 || flag30 || flag31 || flag32 || flag33 || flag34 || flag35 || flag36 || flag37 || flag38 || flag39 || flag40 || flag41 || flag42 || flag43 || flag44 || flag45 || flag46 || flag47 || flag48 || flag49 || flag50 || flag51 || flag52 || flag53 || flag54 || flag55 || flag56 || flag57 || flag58 || flag59 || flag60 || flag61 || flag62 || flag63 || flag64 || flag65 || flag66 || flag67 || flag68 || flag69 || flag71 || flag72 || flag73 || flag74 || flag75 || flag76 || flag77 || flag78 || flag79 || flag80 || flag81 || flag82 || flag83 || flag84 || flag85 || flag86 || flag87 || flag88 || flag89 || flag90 || flag91 || flag92 || flag93 || flag94 || flag95 || flag96){
				 num = num + 1;
                 row[point] = i + 2;
                 col[point] = j + 2;
                 point = point + 1;
			}
			else if(flag70){
				II[(i+2)*256+j+2] = c;
				II[(i+2)*256+j+3] = c;
				II[(i+3)*256+j+2] = c;
				II[(i+3)*256+j+3] = c;
				num = num + 1;
                row[point] = i + 2;
                col[point] = j + 2;
                point = point + 1;
			}
			}
		}
        
        int[] mybin = new int[150*10];
        int[] temp = new int[1000];
        
        for(int p=0;p<point;++p){
        	if(II[row[p]*256+col[p]] == 0){
        		temp[p] = 1;
        	}
        }
        
        Log.i("Point", String.valueOf(point));
        int t = point/10;
        boolean f = true;
        int q = 0;
        for(q=0;q<t;++q){
        	f = false;
        	System.arraycopy(temp,10*q,mybin,10*q,10);
           	for(int i=0;i<10;++i){
           		if(mybin[q*10+i] != 0) {f=true;break;}
           	}
        	
           	if(!f) break;
        }
        
        int[] myinfo = new int[(q-1)*10];
        System.arraycopy(mybin,10,myinfo,0,myinfo.length);
        
        ReedSolomonDecoder decoder = new ReedSolomonDecoder(GenericGF.AZTEC_PARAM);
        int number =  myinfo.length/40;
        int left = (myinfo.length % 40) / 10; 
        
        int[] result = new int[myinfo.length/40*6+left*2];
        int[] tmp = new int[8];
        
        for(int i=0;i<number;++i){
        	tmp[0] = myinfo[40*i+2]*8+myinfo[40*i+3]*4+myinfo[40*i+4]*2+myinfo[40*i+5];
        	tmp[1] = myinfo[40*i+6]*8+myinfo[40*i+7]*4+myinfo[40*i+8]*2+myinfo[40*i+9];
        	tmp[2] = myinfo[40*i+12]*8+myinfo[40*i+13]*4+myinfo[40*i+14]*2+myinfo[40*i+15];
        	tmp[3] = myinfo[40*i+16]*8+myinfo[40*i+17]*4+myinfo[40*i+18]*2+myinfo[40*i+19];
        	tmp[4] = myinfo[40*i+22]*8+myinfo[40*i+23]*4+myinfo[40*i+24]*2+myinfo[40*i+25];
        	tmp[5] = myinfo[40*i+26]*8+myinfo[40*i+27]*4+myinfo[40*i+28]*2+myinfo[40*i+29];
        	tmp[6] = myinfo[40*i+32]*8+myinfo[40*i+33]*4+myinfo[40*i+34]*2+myinfo[40*i+35];
        	tmp[7] = myinfo[40*i+36]*8+myinfo[40*i+37]*4+myinfo[40*i+38]*2+myinfo[40*i+39];
        	try {
				decoder.decode(tmp, 2);
			} catch (ReedSolomonException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
        	System.arraycopy(tmp, 0, result, 6*i, 6);
        }
        
        int h = number*6;
        for(int j = 0;j < left;++j){
        	result[h++] = myinfo[40*number+10*j+2]*8 + myinfo[40*number+10*j+3]*4 + myinfo[40*number+10*j+4]*2 + myinfo[40*number+10*j+5];
        	result[h++] = myinfo[40*number+10*j+6]*8 + myinfo[40*number+10*j+7]*4 + myinfo[40*number+10*j+8]*2 + myinfo[40*number+10*j+9];
        }
        
        char[] Final =  new char[result.length/2];
        for(int i = 0; i < Final.length;++i){
        	Final[i] = (char)(result[2*i]*16+result[2*i+1]);
        }
        Log.i("Final", Arrays.toString(Final));
        Log.i("myinfo", Arrays.toString(myinfo));
        
        return Final;
	}
	
	private int get(int[] pix,int i,int j)
    {
        return pix[imgWidth*i+j];
    }
		
	private boolean erode_area8(int[]pix,int u,int v,int foreground)
    {   if(
           get(pix,u-1,v-1)==foreground&&
           get(pix,u-1,v)==foreground&&
           get(pix,u-1,v+1)==foreground&&
           get(pix,u,v+1)==foreground&&
           get(pix,u,v)==foreground  &&  //中心点
           get(pix,u,v-1)==foreground&&
           get(pix,u+1,v-1)==foreground&&
           get(pix,u+1,v)==foreground&&
           get(pix,u+1,v+1)==foreground  )
        return true;
       else return false;
    }
		
		void set(int[]pix,int i,int j, int value){
			pix[imgWidth*i+j]=value;
        }
		
		public void erode(int[] pix,int []erodePix,int foreground)//8邻域腐蚀
	    {
	        for(int u=0;u<imgHeight-1;u++)//u代表hang
	           for(int v=0;v<imgWidth-1;v++)
	                set(erodePix, u , v , 1-foreground);// 涂满背景

	        for(int u=1;u<imgHeight-2;u++)//u代表hang
	           for(int v=1;v<imgWidth-2;v++)
	               if(erode_area8(pix,u,v,foreground) )
	                   set(erodePix,u,v,foreground);
	               else
	                   set(erodePix,u,v,1-foreground);
	    }
		
		public void dilate(int[] pix,int []dilatePix,int foreground)//8邻域膨胀
	    {
	        for(int u=0;u<imgHeight-1;u++)//u代表hang
	           for(int v=0;v<imgWidth-1;v++)
	               set(dilatePix, u, v , 1-foreground);// 涂满背景

	        for(int u=1;u<imgHeight-2;u++)//u代表hang
	           for(int v=1;v<imgWidth-2;v++)
	           {

	               if (get(pix, u, v) == foreground) {
	                   set(dilatePix, u - 1, v - 1, foreground);
	                   set(dilatePix, u - 1, v, foreground);
	                   set(dilatePix, u - 1, v + 1, foreground);
	                   set(dilatePix, u, v + 1, foreground);
	                   set(dilatePix, u, v, foreground); //中心点
	                   set(dilatePix, u, v - 1, foreground);
	                   set(dilatePix, u + 1, v - 1, foreground);
	                   set(dilatePix, u + 1, v, foreground);
	                   set(dilatePix, u + 1, v + 1, foreground);
	               }
	           }
	    }
	
}
