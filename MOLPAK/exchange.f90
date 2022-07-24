      SUBROUTINE EXCHANGE (X,Y) 
          
      REAL :: X, Y, T
      
      T = X                                                              
      X = Y                                                              
      Y = T                                                              
      RETURN                                                             
      END SUBROUTINE EXCHANGE   
