package edu.sjtu.XiZhang.My_Decoder;


import android.app.Activity;
import android.content.Context;
import android.content.Intent;
import android.content.res.Configuration;
import android.graphics.Bitmap;
import android.graphics.BitmapFactory;
import android.graphics.PixelFormat;
import android.hardware.Camera;
import android.hardware.Camera.PictureCallback;
import android.os.Bundle;
import android.os.Handler;
import android.os.Message;
import android.util.Log;
import android.view.SurfaceHolder;
import android.view.SurfaceView;
import android.view.View;
import android.view.View.OnClickListener;
import android.view.Window;
import android.view.WindowManager;
import android.widget.Button;

public class MainActivity extends Activity implements SurfaceHolder.Callback{
	
	//Define Objects
	private SurfaceView mSurfaceView = null;
	private SurfaceHolder mSurfaceHolder  = null;
	private Camera mCamera = null;
	private boolean storeFlag = false, bIfPreview=false;
	private int mPreviewWidth= 720, mPreviewHeight = 1280;
	private Button capBtn;
	
	//Used to receive decoded data
	private Handler uiHandler;
	private Bundle mBundle;
	private String result = null;
	
	static Context context;
	
	//InitSurfaceView
	private void initSurfaceView(){
		mSurfaceView = (SurfaceView) this.findViewById(R.id.preview_view);
		mSurfaceHolder = mSurfaceView.getHolder();
		mSurfaceHolder.addCallback(MainActivity.this);
		mSurfaceHolder.setFormat(PixelFormat.TRANSPARENT);
	}

	@Override
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		//no title
		requestWindowFeature(Window.FEATURE_NO_TITLE);
		//full screen
		getWindow().setFlags(WindowManager.LayoutParams.FLAG_FULLSCREEN, WindowManager.LayoutParams.FLAG_FULLSCREEN);
		
		setContentView(R.layout.code_decoder_layout);
		
		capBtn = (Button)findViewById(R.id.cap_btn);
		initSurfaceView();
		
		uiHandler = new Handler(){
			@Override
			public void handleMessage(Message msg){
				if(msg.what == 0x123){
					mBundle = msg.getData();
					char[] data = mBundle.getCharArray("rdata");
					result = new String(data);
					
					Intent s = new Intent(MainActivity.this, Result.class);
					Bundle tmp = new Bundle();
					tmp.putString("finalData",result);
					s.putExtras(tmp);
					startActivity(s);
				}
			}
		};
		
        context = getApplicationContext();
		
		capBtn.setOnClickListener(new OnClickListener(){
			@Override
			public void onClick(View v){
				storeFlag = true;
			}
		});	
	}
	
	//Store Image into Disk
	/*public static void writeImageToDisk(byte[] img, String fileName){
		File file = new File(Environment.getExternalStorageDirectory().getPath()+"/DCIM/"+fileName+".PNG");
		if(!file.exists())
			try {
				file.createNewFile();
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		try{
		    FileOutputStream fops = new FileOutputStream(file);
			Bitmap bitmap = BitmapFactory.decodeByteArray(img, 0, img.length);
			bitmap.compress(Bitmap.CompressFormat.PNG, 100, fops);
			fops.flush();
			fops.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}*/
	
	//SurfaceHolder.Callback
	public void surfaceCreated(SurfaceHolder holder){
		mCamera = Camera.open();
		try{
			Log.i("TAG", "SurfaceHolder.Callback：surface Created");
			mCamera.setPreviewDisplay(mSurfaceHolder);
		} catch (Exception ex){
			if(mCamera != null){
				mCamera.release();
				mCamera = null;
			}
			Log.i("TAG initCamera", ex.getMessage());
		}
		
		mCamera.setPreviewCallback(new Camera.PreviewCallback() {
			
			@Override
			public void onPreviewFrame(byte[] data, Camera camera) {
				Log.i("First Data Length", (String.valueOf(data.length)));
				// TODO Auto-generated method stub
				if(storeFlag){
					storeFlag = false;
					mCamera.takePicture(null,null,new PictureCallback(){
						public void onPictureTaken(byte[] data, Camera camera){
							//writeImageToDisk(data,"raw");
							Bitmap bitmap = BitmapFactory.decodeByteArray(data,0,data.length);
							//Log.i("ROW", String.valueOf(bitmap.getWidth()));
							int[] pixels = new int[3264*2448];
							bitmap.getPixels(pixels, 0, 2448, 0, 0, 2448, 3264);
							int[] finalpixels = new int[2448*3264];
							for(int i=0;i<2448;++i){
								for(int j=0;j<3264;++j){
									finalpixels[i*3264+j] = pixels[j*2448+2447-i];
								}
							}
							//Log.i("Second Data Length", (String.valueOf(pixels.length)));
							//Log.i("TAG", "Entering Decode Thread!");
							Decoder_Thread mThread = new Decoder_Thread(finalpixels,uiHandler);
							mThread.start();
						}
					});
					//Log.i("TAG", "Entering Decode Thread!");
					//Decoder_Thread mThread = new Decoder_Thread(data,uiHandler);
					//mThread.start();
					//writeImgToDisk(data,"RawImage.bmp");
				}
			}
		});
	}
	
	public void surfaceChanged(SurfaceHolder holder, int format, int width, int height){
		Log.i("TAG", "SurfaceHolder.Callback：Surface Changed!");
		initCamera();
		mCamera.cancelAutoFocus();
	}
	
	public void surfaceDestroyed(SurfaceHolder holder){
		Log.i("TAG", "SurfaceHolder.Callback：Surface Destroyed");
		if(mCamera != null){
			mCamera.setPreviewCallback(null);
			mCamera.stopPreview();
			bIfPreview = false;
			mCamera.release();
			mCamera = null;
		}
	}
	
	@SuppressWarnings("deprecation")
	private void initCamera(){
		Log.i("TAG", "Initing camera.....");
		
		if (bIfPreview)
    	{
    		mCamera.stopPreview();//stopCamera();
    	}
		
		if(mCamera != null){
			try{
				Camera.Parameters parameters = mCamera.getParameters();
				parameters.setPictureFormat(PixelFormat.JPEG);
				parameters.setPreviewFormat(PixelFormat.YCbCr_420_SP);
				/*List<Size> pictureSizes = mCamera.getParameters().getSupportedPictureSizes();
		        List<Size> previewSizes = mCamera.getParameters().getSupportedPreviewSizes();
		        List<Integer> previewFormats = mCamera.getParameters().getSupportedPreviewFormats();
		        List<Integer> previewFrameRates = mCamera.getParameters().getSupportedPreviewFrameRates();
		        Log.i("TAG"+"initCamera", "cyy support parameters is ");
		        Size psize = null;
		        for (int i = 0; i < pictureSizes.size(); i++)
		        {
		        	psize = pictureSizes.get(i);
		        	Log.i("TAG"+"initCamera", "PictrueSize,width: " + psize.width + " height" + psize.height);
		        }
		        for (int i = 0; i < previewSizes.size(); i++)
		        {
		        	psize = previewSizes.get(i);
		        	Log.i("TAG"+"initCamera", "PreviewSize,width: " + psize.width + " height" + psize.height);
		        }
		        Integer pf = null;
		        for (int i = 0; i < previewFormats.size(); i++)
		        {
		        	pf = previewFormats.get(i);
		        	Log.i("TAG"+"initCamera", "previewformates:" + pf);
		        }*/

		        /*List<int[]> supportedPreviewFpsRange = mCamera.getParameters().getSupportedPreviewFpsRange();
		        List<Integer> supportedPreviewFrameRates = mCamera.getParameters().getSupportedPreviewFrameRates();
		        int fr;
		        for (int i = 0; i < supportedPreviewFrameRates.size(); i++)
			    {
		        	fr = supportedPreviewFrameRates.get(i);
			       	Log.i("TAG"+"initCamera", "previewformates:" + fr);
			    }*/
				parameters.setFocusMode(Camera.Parameters.FOCUS_MODE_CONTINUOUS_PICTURE);
				parameters.setPictureSize(3264, 2448);
				parameters.setPreviewSize(mPreviewHeight, mPreviewWidth);
				
				// 横竖屏镜头自动调整
		        if (this.getResources().getConfiguration().orientation != Configuration.ORIENTATION_LANDSCAPE) 
		        {
		        	parameters.set("orientation", "portrait"); //
		        	parameters.set("rotation", 90); // 镜头角度转90度（默认摄像头是横拍） 
		        	mCamera.setDisplayOrientation(90); // 在2.2以上可以使用
		        } else// 如果是横屏
		        {
		        	parameters.set("orientation", "landscape"); //
		        	mCamera.setDisplayOrientation(0); // 在2.2以上可以使用
		        } 
		        
		        // 设定配置参数并开启预览
		        mCamera.setParameters(parameters); // 将Camera.Parameters设定予Camera    
		        mCamera.startPreview(); // 打开预览画面
		        mCamera.cancelAutoFocus();// 2如果要实现连续的自动对焦，这一句必须加上 
		        bIfPreview = true;
		        
		        // 【调试】设置后的图片大小和预览大小以及帧率
		        /*Camera.Size csize = mCamera.getParameters().getPreviewSize();
		        mPreviewHeight = csize.height; //
		        mPreviewWidth = csize.width;
		        Log.i("TAG"+"initCamera", "after setting, previewSize:width: " + csize.width + " height: " + csize.height);
		        csize = mCamera.getParameters().getPictureSize();
		        Log.i("TAG"+"initCamera", "after setting, pictruesize:width: " + csize.width + " height: " + csize.height);
		        Log.i("TAG"+"initCamera", "after setting, previewformate is " + mCamera.getParameters().getPreviewFormat());
		        Log.i("TAG"+"initCamera", "after setting, previewframerate is " + mCamera.getParameters().getPreviewFrameRate());*/
			} catch (Exception e)
	    	{ 
		    	   e.printStackTrace();
		    	}
		}
	}

	
}
