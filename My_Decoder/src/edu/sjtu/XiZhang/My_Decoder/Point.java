package edu.sjtu.XiZhang.My_Decoder;


public class Point {		//一个point类,包含了横纵坐标
	
	public int x;
	public int y;
	
	Point()
	{
		x=0;
		y=0;
	}
	
	Point(int x, int y)
	{
		this.x = x;
		this.y = y;
	}
	
	public void set(int xi, int yi)
	{
		x = xi;
		y = yi;
	}

	public Point middlePoint(Point point0){
		int tmpx=(x+point0.x)/2;
		int tmpy=(y+point0.y)/2;
		return new Point(tmpx,tmpy);
	}
}