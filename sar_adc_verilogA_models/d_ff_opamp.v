//Amp_Diff, veriloga
//2019.11.04  by Hu Bingxiang

`include "constants.vams"
`include "disciplines.vams"
`define PI 3.1415926535

module diff_opamp(voutp,voutn,vref,vinp,vinn,vsupplyp,vsupplyn);
input vref, vsupplyp, vsupplyn;
inout voutp, voutn, vinp, vinn;

parameter real gain = 835e1;
parameter real freq_unitygain = 1.0e8;
parameter real rin = 1e6;
parameter real vin_offset = 0.0;
parameter real ibias = 0.0;
parameter real iin_max = 1e-6;
parameter real slew_rate = 0.5e6;
parameter real rout = 80;
parameter real vsoft = 0.5;
parameter real vmax_in = 2;
	real c1;
	real gm_nom; 
	real r1;
 	real vin_val;
 
electrical voutp, voutn, vref, vinp, vinn, vsupplyp, vsupplyn;
electrical coutp, coutn;
	analog 
	begin
		@(initial_step or initial_step("dc", "ac", "tran", "xf")) 
		begin
			c1 = iin_max / (slew_rate);
			gm_nom = 2 * `PI * freq_unitygain * c1;
			r1 = gain / gm_nom;
	//		vmax_in = 2;
		end
	
		vin_val = V(vinp,vinn)/2 + vin_offset;
		//input stage
		I(vinp,vinn) <+ (V(vinp,vinn) + vin_offset)/rin;
		I(vref,vinp) <+ ibias;
		I(vref,vinn) <+ ibias;
		
		//GM stage with slewing
		I(vref,coutp) <+ V(vref,coutp)/100e6;
		I(vref,coutn) <+ V(vref,coutn)/100e6;
		
		if(vin_val > vmax_in) 
			begin
				I(vref,coutp) <+ iin_max;
				I(vref,coutn) <+ iin_max;
			end
		
		else if (vin_val < -vmax_in) 
			begin
				I(vref,coutp) <+ -iin_max;
				I(vref,coutn) <+ iin_max;
			end
		else begin
			I(vref,coutp) <+ 0.5 * gm_nom * vin_val;
			I(vref,coutn) <+ -0.5 * gm_nom * vin_val;
			end
		
		//Dominant pole
		I(coutp,vref) <+ ddt(c1 * V(coutp, vref));
		I(coutp,vref) <+ V(coutp, vref)/r1;
		I(coutn,vref) <+ ddt(c1 * V(coutn, vref));
		I(coutn,vref) <+ V(coutn, vref)/r1;
		
		//Output Stage
		I(vref, voutp)<+ V(coutp, vref) / rout;
		I(voutp, vref)<+ V(voutp, vref) / rout;
		I(vref, voutn)<+ V(coutn, vref) / rout;
		I(voutn, vref)<+ V(voutn, vref) / rout;

		//Soft Output Limiting
		if  ( V(voutp) > (V(vsupplyp)- vsoft) )
			I(coutp, vref)<+ gm_nom *(V(voutp, vsupplyp)+ vsoft);
		else if ( V(voutp) < (V(vsupplyn)- vsoft) )
			I(coutp, vref)<+ gm_nom *(V(voutp, vsupplyn)- vsoft);
		
		if ( V(voutn) > (V(vsupplyp)- vsoft) )
			I(coutn, vref) <+ gm_nom *(V(voutn, vsupplyp)+ vsoft);
		else if ( V(voutn) < (V(vsupplyn)- vsoft) )
			I(coutn, vref) <+ gm_nom *(V(voutn, vsupplyn)- vsoft);
	end

endmodule
