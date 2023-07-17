TH1D *toGorH(TH1D *hist);
TH1D *toGorH_inverse(TH1D *hist);

void draw_test(TCanvas *cvs, int seed, int iseed)
{
    bool draw_process1 = false;
    bool draw_process2 = false;

    bool test_setting = false;
    bool paper_setting = true;

    int    num_deconvolution = 20;
    int    num_events = 10000; // 10000
    int    num_particles = 10; // 10
    int    num_bins = 101;
    double x_max =  2.;
    double x_min = -2.;
    double x_range = 1;
    double x_range_off = 0;
    double lambda = 0.01;
    double acc_power = 1.99;

    double hist_min = 0;
    double hist_max = 0;
    double factor_height = 1;
    int width = 3;

    if (test_setting) {
        num_deconvolution = 1;
        num_events = 10000;
        num_particles = 10;
        num_bins = 8;
        x_max =  4.;
        x_min = -4.;
        x_range = 1;
        lambda = 0.01;
        acc_power = 1.99;
        hist_min = x_min;
        hist_max = x_max;
        factor_height = 1.5;
        width = 5;
    }
    else if (paper_setting) {
        num_deconvolution = 2200;
        num_events = 10000;
        num_particles = 10;
        num_bins = 41;
        x_max =  2.05;
        x_min = -2.05;
        x_range = 1;
        lambda = 0.01;
        acc_power = 1.99;
        width = 5;
    }

    double bin_width = (x_max-x_min)/num_bins;
    if (hist_min==0) hist_min = x_min * 0.75;
    if (hist_max==0) hist_max = x_max * 0.75;
    if (hist_min > -x_range*1.2) hist_min = 1.2*(-x_range);
    if (hist_max < +x_range*1.2) hist_max = 1.2*(+x_range);
    double height = factor_height * 2 * num_particles / (x_range*2.);

    int num_off = 0;
    if (num_bins%2==0)
        num_off = num_bins/2;
    else
        num_off = (num_bins+1)/2;

    auto hist_x_measured        = new TH1D(Form("hist_x_measured_%d",iseed)     ,"x measured",                num_bins,x_min,x_max);
    auto hist_x_init            = new TH1D(Form("hist_x_init_%d",iseed)         ,"x init /bin/event",         num_bins,x_min,x_max);
    auto hist_x_smearing        = new TH1D(Form("hist_x_smearing_%d",iseed)     ,"x smearing",                num_bins,x_min,x_max);
    auto hist_x_in_process_2    = new TH1D(Form("hist_x_in_process_2_%d",iseed) ,"convolution of smearing and measured", num_bins,x_min,x_max);
    auto hist_x_in_process_3    = new TH1D(Form("hist_x_in_process_3_%d",iseed) ,"x in deblurring process 3", num_bins,x_min,x_max);
    auto hist_x_restored        = new TH1D(Form("hist_x_restored_%d",iseed)     ,"x restored",                num_bins,x_min,x_max);

    double x_mean_average_error_squared = 0.;
    for (auto i_event=0; i_event<num_events; ++i_event) {
        double x_mean=0.;
        double x_squared_mean=0.;
        double x_init_array[num_particles]; 
        for (auto i_particle=0; i_particle<num_particles; ++i_particle)
        {
            //double x_init = gRandom -> Uniform(-x_range, x_range);
            double x_init = gRandom -> Uniform(-x_range+x_range_off, x_range+x_range_off);
            x_init_array[i_particle] = x_init;
            x_mean         += x_init;
            x_squared_mean += x_init*x_init;
        }
        x_mean            = x_mean         / num_particles;
        x_squared_mean    = x_squared_mean / num_particles;
        //!estimate of average error squared
        double x_variance = (x_squared_mean-x_mean*x_mean);
        double x_average_error_squared = (x_squared_mean-x_mean*x_mean) / (num_particles-1);
        {
            double x_variance2 = 0;
            for (auto i_particle=0; i_particle<num_particles; ++i_particle) {
                x_variance2 += (x_mean-x_init_array[i_particle])*(x_mean-x_init_array[i_particle]);
            }
            x_variance2 = x_variance2/(num_particles-1);
            //cout << x_average_error_squared << " " << x_variance2 << endl;
        }
        x_mean_average_error_squared += x_average_error_squared;
        for (auto i_particle=0; i_particle<num_particles; ++i_particle)
        {
            double x_init = x_init_array[i_particle];
            hist_x_init -> Fill(x_init);
            double x_normalized = (num_particles*x_mean - x_init) / (num_particles-1);
            double x_renormalized = x_init - x_normalized;
            hist_x_measured -> Fill(x_renormalized);
            hist_x_restored -> Fill(x_renormalized);
        }
    }

    /// renormalized width of average position: ptcle vs rest
    x_mean_average_error_squared = x_mean_average_error_squared / num_events;

    //SIGXK=XAKS*NPART/(NPART-1)  !renormalized width of average position: ptcle vs rest
    double sigmaX = x_mean_average_error_squared * num_particles/(num_particles-1);
    double sqrt_vn = TMath::Sqrt(sigmaX);
    double norm_factor1 = 1. / (sqrt_vn * TMath::Sqrt(2*TMath::Pi()));
    cout << endl;
    cout << x_mean_average_error_squared << " " << sigmaX << " " << sqrt_vn << " " << norm_factor1 << endl;

    {
        double a_average_error_squared = (x_range*2)*(x_range*2)/12.; // variance
        double a_mean_average_error_squared = a_average_error_squared / num_particles;

        double sigmaA = a_mean_average_error_squared * num_particles / (num_particles-1);
        double sqrt_vn2 = TMath::Sqrt(sigmaA);
        double norm_factor2 = 1. / (sqrt_vn2 * TMath::Sqrt(2*TMath::Pi()));
        cout << a_mean_average_error_squared << " " << sigmaA << " " << sqrt_vn2 << " " << norm_factor2 << endl;
        cout << endl;

        sigmaX = sigmaA;
        norm_factor1 = norm_factor2;
    }


    cout << "num_deconvolution : " << num_deconvolution << endl;
    cout << "num_events          " << num_events << endl;
    cout << "num_particles       " << num_particles << endl;
    cout << "num_bins            " << num_bins << endl;
    cout << "x_max               " << x_max << endl;
    cout << "x_min               " << x_min << endl;
    cout << "x_range             " << x_range << endl;
    cout << "lambda              " << lambda << endl;
    cout << "acc_power           " << acc_power << endl;
    cout << "hist_min            " << hist_min << endl;
    cout << "hist_max            " << hist_max << endl;
    cout << "factor_height       " << factor_height << endl;
    cout << "width               " << width << endl;
    cout << endl;

    //cout << "variance_uniform_dist      : " << variance_uniform_dist  << endl;
    //cout << "variance_div_num_particles : " << variance_div_num_particles  << endl;
    //cout << "sigma_new                  : " << sigma_new  << endl;
    //cout << "norm_factor                : " << norm_factor  << endl;
    cout << "x_mean_average_error_squared         : " << x_mean_average_error_squared  << endl;
    cout << "sigmaX                     : " << sigmaX  << endl;
    //cout << "sqrt_vn                    : " << sqrt_vn  << endl;
    cout << "norm_factor1               : " << norm_factor1  << endl;
    cout << endl;

    for (auto hbin=1; hbin<=num_bins; ++hbin) {
        double bin_x = hist_x_smearing -> GetBinCenter(hbin);
        double smearing = norm_factor1 * TMath::Exp(-.5*bin_x*bin_x/sigmaX)*bin_width;
        hist_x_smearing -> SetBinContent(hbin,smearing);
    }

    double sum_x_renormalized=0.;
    for (auto hbin=1; hbin<=num_bins; ++hbin)
        sum_x_renormalized += hist_x_measured -> GetBinContent(hbin);

    int i_pad_restore = 0;
    TCanvas *cvs_restore = nullptr;
    if (draw_process1) {
        cvs_restore = new TCanvas(Form("cvs_restore_%d",iseed),"",1200,750);
        cvs_restore -> Divide(4,3);
        //for (auto i=1; i<=12; ++i) cvs_restore -> cd(i) -> SetGrid();
    }

    /// deconvolution
    for (auto i_iterate=1; i_iterate<=num_deconvolution; ++i_iterate)
    {
        for (auto hbin2=1; hbin2<=num_bins; ++hbin2)
        {
            int count_added = 0;
            double x_smear_restore = 0.;
            for (auto hbin_s=1; hbin_s<=num_bins; ++hbin_s)
            {
                int hbin_r = hbin2 - hbin_s + num_off;
                if (hbin_r>0 && hbin_r<=num_bins) {
                    count_added++;
                    double x_smeared = hist_x_smearing -> GetBinContent(hbin_s);
                    double x_restore = hist_x_restored -> GetBinContent(hbin_r);
                    x_smear_restore += x_smeared * x_restore;
                }
            }
            /// convoluted current restore w/blurring function
            hist_x_in_process_2 -> SetBinContent(hbin2,x_smear_restore);
        }
        //break;

        for (auto hbin_r=1; hbin_r<=num_bins; ++hbin_r)
        {
            double x_restore1=0.;
            for (auto hbin_s=1; hbin_s<=num_bins; ++hbin_s)
            {
                int hbin_m = hbin_r + hbin_s - num_off;
                if (hbin_m>0 && hbin_m<=num_bins)
                {
                    double x_smeared        = hist_x_smearing -> GetBinContent(hbin_s);
                    double x_renormalized   = hist_x_measured -> GetBinContent(hbin_m);
                    double x_smear_restore  = hist_x_in_process_2 -> GetBinContent(hbin_m);
                    if (x_smear_restore>0.)
                        x_restore1 += x_smeared * x_renormalized / x_smear_restore;
                    else
                        x_restore1 += x_smeared;
                }
            }


            double regularization_factor = 1.;
            if (hbin_r>1 && hbin_r<=num_bins-1)
            {
                double x_reno_0 = hist_x_restored -> GetBinContent(hbin_r);
                double x_reno_b = hist_x_restored -> GetBinContent(hbin_r-1);
                double x_reno_a = hist_x_restored -> GetBinContent(hbin_r+1);
                if      (x_reno_0 > x_reno_b && x_reno_0 > x_reno_a) regularization_factor = 1./(1. + lambda); // peak
                else if (x_reno_0 < x_reno_b && x_reno_0 < x_reno_a) regularization_factor = 1./(1. - lambda); // bump
            }

            double x_restore2 = hist_x_restored -> GetBinContent(hbin_r);
            double x_restore3 = x_restore2 * regularization_factor * TMath::Power(x_restore1, acc_power);
            hist_x_in_process_3 -> SetBinContent(hbin_r, x_restore3);
        }

        double sum_of_xsr=0.;
        for (auto hbin=1; hbin<=num_bins; ++hbin)
            sum_of_xsr += hist_x_in_process_3 -> GetBinContent(hbin);

        for (auto hbin=1; hbin<=num_bins; ++hbin) {
            double content = hist_x_in_process_3 -> GetBinContent(hbin) * sum_x_renormalized / sum_of_xsr;
            hist_x_restored -> SetBinContent(hbin,content);
        }

        if (draw_process1) {
            if (
                    (num_deconvolution<=15&&i_iterate<15) ||
                    (i_iterate==1 ||
                     i_iterate==2 ||
                     i_iterate==3 ||
                     i_iterate==5 ||
                     i_iterate==10 ||
                     i_iterate==50 ||
                     i_iterate==100 ||
                     i_iterate==200 ||
                     i_iterate==500 ||
                     i_iterate==1000 ||
                     i_iterate==1500 ||
                     i_iterate==2000 ||
                     i_iterate==3000 ||
                     i_iterate==4000 ||
                     i_iterate==5000)
               )
            {
                ++i_pad_restore;
                cvs_restore -> cd(i_pad_restore);
                auto hist_restored_clone = (TH1D *) hist_x_restored -> Clone(Form("restored_it%d_%d",i_iterate,iseed));
                hist_restored_clone -> SetTitle(Form("restored with iteration=%d",i_iterate));
                hist_restored_clone -> Draw();
            }
        }
    }

    hist_x_init -> SetMinimum(0);
    hist_x_init -> Scale(1./(bin_width*num_events));
    hist_x_measured -> Scale(1./(bin_width*num_events));
    hist_x_restored -> Scale(1./(bin_width*num_events));
    hist_x_in_process_2 -> Scale(1./(bin_width*num_events));
    hist_x_in_process_3 -> Scale(1./(bin_width*num_events));

    hist_x_in_process_2 -> SetLineColor(kPink-9);
    hist_x_in_process_2 -> SetLineWidth(3);

    hist_x_measured -> SetLineColor(kBlue-4);
    hist_x_measured -> SetLineWidth(3);

    hist_x_smearing -> SetLineColor(kBlue-4);
    hist_x_smearing -> SetLineWidth(3);

    hist_x_restored -> SetLineColor(kBlack);
    hist_x_restored -> SetLineWidth(3);

    if (draw_process2) {
        auto cvs_all = new TCanvas(Form("cvs_all_%d",iseed),"",1200,700);
        //cvs_all -> Divide(3,2);
        cvs_all -> Divide(2,2);
        cvs_all -> cd(1); hist_x_measured     -> Draw("hist");
        //cvs_all -> cd(3); hist_x_init         -> Draw("hist");
        cvs_all -> cd(2); hist_x_smearing     -> Draw("hist");
        cvs_all -> cd(3); hist_x_in_process_2 -> Draw("hist");
        //cvs_all -> cd(5); hist_x_in_process_3 -> Draw("hist");
        cvs_all -> cd(4); hist_x_restored     -> Draw("hist");
        cvs_all -> SaveAs("figures/summary.png");
    }

    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;
    double sum5 = 0;
    double sum6 = 0;
    int i1 = (1+num_bins/8*2);
    int i2 = (1+num_bins/8*3);
    int i3 = (num_bins-num_bins/8*3);
    int i4 = (num_bins-num_bins/8*2);
    double x1 = hist_x_measured -> GetXaxis() -> GetBinLowEdge(i1);
    double x2 = hist_x_measured -> GetXaxis() -> GetBinUpEdge(i2);
    double x3 = hist_x_measured -> GetXaxis() -> GetBinLowEdge(i3);
    double x4 = hist_x_measured -> GetXaxis() -> GetBinUpEdge(i4);
    cout << "i: " << i2 << " -> " << i1 << endl;
    cout << "x: " << x2 << " -> " << x1 << endl;
    cout << "i: " << i3 << " -> " << i4 << endl;
    cout << "x: " << x3 << " -> " << x4 << endl;
    for (auto i=i1; i<=i2; ++i) {
        sum1 += hist_x_measured -> GetBinContent(i);
        sum2 += hist_x_init -> GetBinContent(i);
        sum3 += hist_x_restored -> GetBinContent(i);
    }
    for (auto i=i3; i<=i4; ++i) {
        sum4 += hist_x_measured -> GetBinContent(i);
        sum5 += hist_x_init -> GetBinContent(i);
        sum6 += hist_x_restored -> GetBinContent(i);
    }
    cout << "measured: " << sum1 << " " << sum4 << endl;
    cout << "original: " << sum2 << " " << sum5 << endl;
    cout << "restored: " << sum3 << " " << sum6 << endl;

    TString draw_option = "histsame";

    auto hist1 = toGorH(hist_x_init);
    auto hist2 = toGorH(hist_x_measured);
    auto hist5 = toGorH_inverse(hist_x_measured);
    auto hist3 = toGorH(hist_x_restored);
    auto hist4 = toGorH(hist_x_in_process_2);

    cvs -> cd(iseed);
    auto hist = new TH2D(Form("frame_%d",iseed),";x;dN/dx;",100,hist_min,hist_max,100,0,height);
    hist -> SetStats(0);
    hist -> Draw();

    hist1 -> SetLineColor(kPink-9);
    hist1 -> SetLineWidth(width);
    hist1 -> DrawClone(draw_option);

    hist2 -> SetLineColor(kBlue-4);
    hist2 -> SetLineWidth(width);
    hist2 -> SetLineStyle(9);
    hist2 -> DrawClone(draw_option);

    hist5 -> SetLineColor(kMagenta-9);
    hist5 -> SetLineWidth(width);
    hist5 -> SetLineStyle(1);
    //hist -> DrawClone(draw_option);

    hist3 -> SetLineColor(kGreen+1);
    hist3 -> SetLineWidth(width);
    hist3 -> SetLineStyle(1);
    hist3 -> DrawClone(draw_option);

    hist4 -> SetLineColor(kCyan);
    hist4 -> SetLineWidth(width);
    hist4 -> SetLineStyle(9);
    //hist4 -> DrawClone(draw_option);

    auto legend = new TLegend();
    legend -> AddEntry(hist1,"original");
    legend -> AddEntry(hist2,"measured");
    legend -> AddEntry(hist3,"restored");
    legend -> SetX1(.14);
    legend -> SetX2(.48);
    legend -> SetY1(.60);
    legend -> SetY2(.88);
    legend -> SetTextSize(0.038);
    legend -> SetBorderSize(1);
    legend -> Draw();

    auto legend2 = new TLegend();
    auto marker = new TMarker(0,0,20);
    legend2 -> AddEntry(marker,Form("# events = %d",num_events),"p");
    legend2 -> AddEntry(marker,Form("binning = (%d, %.2f, %.2f)",num_bins,x_min,x_max),"p");
    legend2 -> AddEntry(marker,Form("# iteration = %d",num_deconvolution),"p");
    legend2 -> AddEntry(marker,Form("# particles = %d",num_particles),"p");
    legend2 -> AddEntry(marker,Form("acc. power = %.2f",acc_power),"p");
    legend2 -> AddEntry(marker,Form("#lambda = %.2f",lambda),"p");
    legend2 -> SetX1(.52);
    legend2 -> SetX2(.86);
    legend2 -> SetY1(.60);
    legend2 -> SetY2(.88);
    legend2 -> SetTextSize(0.038);
    legend2 -> SetMargin(0.15);
    legend2 -> SetBorderSize(1);
    legend2 -> Draw();
}

TH1D *toGorH(TH1D *hist) { return hist; }
TH1D *toGorH_inverse(TH1D *hist)
{
    auto hist_inverse = (TH1D *) hist -> Clone();
    auto nbins = hist -> GetXaxis() -> GetNbins();
    for (auto bin=1; bin<=nbins; ++bin) {
        auto bin2 = (nbins-bin)+1;
        hist_inverse -> SetBinContent(bin2,hist->GetBinContent(bin));
    }
    return hist_inverse;
}

void test_1D() {
    int seed = time(0);
    if (1) {
        auto cvs = new TCanvas();
        cvs -> SetGridy();
        //cvs -> SetGridx();
        gRandom -> SetSeed(seed);
        cout << "seed: " << seed << endl;
        draw_test(cvs,seed,0);
        cvs -> SaveAs("figures/test.png");
        //cvs -> SaveAs("figures/test.eps");
        //cvs -> SaveAs("figures/test.pdf");
    }
    else {
        auto cvs = new TCanvas("cvs","",3600,2000);
        cvs -> Divide(3,2);
        for (auto iseed=1; iseed<=6; ++iseed) {
            seed = seed + iseed;
            gRandom -> SetSeed(seed);
            cout << "seed: " << seed << endl;
            draw_test(cvs,seed,iseed);
        }
    }
}
