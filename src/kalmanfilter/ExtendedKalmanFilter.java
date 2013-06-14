package kalmanfilter;

/**
 *
 * @author nick
 */
public class ExtendedKalmanFilter {
    
    Matrix x;
    
    public ExtendedKalmanFilter(double initialX, double initialY, double initialHeading) {
        x = new Matrix(new double[][] { {initialX}, {initialY}, {initialHeading},
            {0}, {0}, {0}});
    }

    Matrix u = new Matrix(new double[][] { {0}, {0}, {0}, {0}, {0}, {0} }); 
    
    Matrix P = new Matrix(new double[][] { {0, 0, 0, 0, 0, 0},
                                           {0, 0, 0, 0, 0, 0},
                                           {0, 0, 0, 0, 0, 0},
                                           {0, 0, 0, 1000, 0, 0},
                                           {0, 0, 0, 0, 1000, 0},
                                           {0, 0, 0, 0, 0, 1000} });
    
    Matrix K = new Matrix(new double[][] { {0, 0},
                                           {0, 0},
                                           {0, 0},
                                           {0, 0},
                                           {0, 0},
                                           {0, 0} });
 
    Matrix H = new Matrix(new double[][] { {1, 0, 0, 0, 0, 0},
                                           {0, 1, 0, 0, 0, 0},
                                           {0, 0, 1, 0, 0, 0} });
    
    Matrix R = new Matrix(new double[][] { {.1, 0, 0},
                                           {0, .1, 0}, 
                                           {0, 0, .1} });
    
    Matrix I = new Matrix(new double[][] { {1, 0, 0, 0, 0, 0},
                                           {0, 1, 0, 0, 0, 0},
                                           {0, 0, 1, 0, 0, 0},
                                           {0, 0, 0, 1, 0, 0},
                                           {0, 0, 0, 0, 1, 0},
                                           {0, 0, 0, 0, 0, 1} });
    
    public Matrix[] filter(double measure_x, double measure_y, double measure_heading, double dt) {
        
        Matrix F = new Matrix(new double[][] { {1, 0, 0, dt, 0, 0},
                                               {0, 1, 0, 0, dt, 0},
                                               {0, 0, 1, 0, 0, dt},
                                               {0, 0, 0, 1, 0, 0},
                                               {0, 0, 0, 0, 1, 0},
                                               {0, 0, 0, 0, 0, 1} });
        
        x = F.times(this.x).plus(this.u);
        P = (F.times(P).times(F.transpose())).plus(this.u); // ...a little messed up
        
        Matrix Z = new Matrix(new double[][] { { measure_x, measure_y, measure_heading } });
        
        Matrix y = Z.transpose().minus(H.times(x));
        
        Matrix S = H.times(P).times(H.transpose()).plus(R); 
        
        Matrix K = P.times(H.transpose()).times(S.inverse());
        
        x = x.plus(K.times(y));
        this.P = (this.I.minus((K.times(H)))).times(this.P);
        
        return new Matrix[] { x, P };
    }
    
    public static void main(String args[]) {
        
        ExtendedKalmanFilter filter = new ExtendedKalmanFilter(1.0, 0.0, 0.0);
        
        double[][] measurements = { {2, 30, 30}, {3, 60, 60}, {4, 90, 90}, {5, 120, 120},
                                    {6, 150, 150}, {7, 180, 180}};
        
        for (int i = 0; i < measurements.length; i++) {
            filter.filter(measurements[i][0], measurements[i][1], measurements[i][2], .1);
        }
        
        filter.x.print(0, 2);
        
    }
    
}