package edu.sjtu.XiZhang.My_Decoder;

import android.content.Context;
import android.content.res.Resources;
import android.graphics.Bitmap;
import android.graphics.Canvas;
import android.graphics.Color;
import android.graphics.Paint;
import android.graphics.Rect;
import android.graphics.Typeface;
import android.util.AttributeSet;
import android.view.View;

public final class ViewFinderView extends View{
	
	private static final int OPAQUE = 0xFF;
	private static final int CORNER_WIDTH = 10;
	private static final int TEXT_SIZE = 16;
	private static final int TEXT_PADDING_TOP = 100;
	
	private static float density;
	private int screenRate;
	private Paint paint;
	private final int maskColor;
	private final int resultColor;
	private Bitmap resultBitmap;
	
	public ViewFinderView(Context context, AttributeSet attrs){
		super(context, attrs);
		
		//get the density of phone screen
		density = context.getResources().getDisplayMetrics().density;
		
		//transfer pixels to dps
		screenRate = (int)(20*density);
		
		paint = new Paint();
		Resources resources = getResources();
		maskColor = resources.getColor(R.color.viewfinder_mask);
		resultColor =  resources.getColor(R.color.result_view);
		
	}
	
	@Override
	public void onDraw(Canvas canvas){
		int width =  canvas.getWidth();
		int height = canvas.getHeight();
		Rect frame = new Rect(width/2-250,height/2-250,width/2+250,height/2+250);
		if(frame == null) return;
		paint.setColor(resultBitmap != null ? resultColor : maskColor);
		
		//Draw the areas outside of the Rect
		canvas.drawRect(0,0,width,frame.top,paint);
		canvas.drawRect(0, frame.top, frame.left, frame.bottom + 1, paint);  
        canvas.drawRect(frame.right + 1, frame.top, width, frame.bottom + 1,  
                paint);  
        canvas.drawRect(0, frame.bottom + 1, width, height, paint);  
        if (resultBitmap != null) {  
            // Draw the opaque result bitmap over the scanning rectangle  
            paint.setAlpha(OPAQUE);  
            canvas.drawBitmap(resultBitmap, frame.left, frame.top, paint);  
        } else {  
  
            //画扫描框边上的角，总共8个部分  
            paint.setColor(Color.GREEN);  
            canvas.drawRect(frame.left, frame.top, frame.left + screenRate,  
                    frame.top + CORNER_WIDTH, paint);  
            canvas.drawRect(frame.left, frame.top, frame.left + CORNER_WIDTH, frame.top  
                    + screenRate, paint);  
            canvas.drawRect(frame.right - screenRate, frame.top, frame.right,  
                    frame.top + CORNER_WIDTH, paint);  
            canvas.drawRect(frame.right - CORNER_WIDTH, frame.top, frame.right, frame.top  
                    + screenRate, paint);  
            canvas.drawRect(frame.left, frame.bottom - CORNER_WIDTH, frame.left  
                    + screenRate, frame.bottom, paint);  
            canvas.drawRect(frame.left, frame.bottom - screenRate,  
                    frame.left + CORNER_WIDTH, frame.bottom, paint);  
            canvas.drawRect(frame.right - screenRate, frame.bottom - CORNER_WIDTH,  
                    frame.right, frame.bottom, paint);  
            canvas.drawRect(frame.right - CORNER_WIDTH, frame.bottom - screenRate,  
                    frame.right, frame.bottom, paint);
            
            //画扫描框下面的字  
            paint.setColor(Color.WHITE);  
            paint.setTextSize(TEXT_SIZE * density);  
            paint.setAlpha(0x40);  
            paint.setTypeface(Typeface.create("System", Typeface.BOLD));  
            canvas.drawText(getResources().getString(R.string.scan_text), frame.left-100, (float) (frame.bottom + (float)TEXT_PADDING_TOP), paint);
        }
	}

}
