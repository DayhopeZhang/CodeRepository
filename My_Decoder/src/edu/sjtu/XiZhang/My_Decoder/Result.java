package edu.sjtu.XiZhang.My_Decoder;

import android.app.Activity;
import android.os.Bundle;
import android.widget.TextView;

public class Result extends Activity {
	Bundle mBundle;
	String result;
	TextView text_view;
	
	@Override
	protected void onCreate(Bundle savedInstanceState){
		super.onCreate(savedInstanceState);
		setContentView(R.layout.activity_main);
		text_view = (TextView) this.findViewById(R.id.text_view);
		mBundle = getIntent().getExtras();
		result = mBundle.getString("finalData");
		text_view.setText(result);
		
		//setContentView(R.layout.activity_main);
	}
	

}
