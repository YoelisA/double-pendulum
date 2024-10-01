module Main where

import Graphics.Gloss
import Graphics.Gloss.Interface.Pure.Simulate

-- Pendulum parameters
l1, l2 :: Float
l1 = 30    -- Length of the first rod
l2 = 30    -- Length of the second rod

m1, m2 :: Float
m1 = 5     -- Mass of the first bob (used for rendering size)
m2 = 5     -- Mass of the second bob (used for rendering size)


-- Main function
main :: IO ()
main = simulate
  -- Display an empty window
    (InWindow "Double drawDoublePendulum" (400, 400) (10, 10)) -- Window title, size, and position
    white
    60
    initialState                                           -- Background color
    drawPendulum
    updatePendulum                                           -- Picture to display (Blank means nothing)

-- Initial state (angles and angular velocities)
initialState :: (Float, Float, Float, Float)
initialState = (pi / 4, pi / 2, 0, 0)  -- (angle1, angle2, velocity1, velocity2)


drawPendulum :: (Float, Float, Float, Float) -> Picture
drawPendulum (t1, t2, o1, o2) = pictures
  [ color black $ circleSolid 5
  , line [(0, 0), bob1Position]                -- First rod
  , line [bob1Position, bob2Position]          -- Second rod
  , translate (fst bob1Position) (snd bob1Position) (color red $ circleSolid m1) -- First bob
  , translate (fst bob2Position) (snd bob2Position) (color green $ circleSolid m2) -- Second bob
  -- , textDisplay
  ] where
    bob1Position = (l1 * sin t1, - (l1 * cos t1))
    bob2Position = let (x1, y1) = bob1Position
                   in (x1 + l2 * sin t2, y1 - l2 * cos t2)
    textDisplay = Translate (-3) 2 $ Scale 0.1 0.1 $
      Color black $ Text $ 
        "Theta1: " ++ show t1 ++ "\n" ++
        "Theta2: " ++ show t2 ++ "\n" ++
        "Omega1: " ++ show o1 ++ "\n" ++
        "Omega2: " ++ show o2 ++ "\n"
g :: Float
g = 9.81

-- Update function using the equations of motion for a double pendulum
updatePendulum :: ViewPort -> Float -> (Float, Float, Float, Float) -> (Float, Float, Float, Float)
updatePendulum _ dt state =
  rk4Step dt state

-- Function to calculate the derivatives (right-hand side of the ODEs)
derivatives' :: (Float, Float, Float, Float) -> (Float, Float, Float, Float)
derivatives' (theta1, theta2, omega1, omega2) =
  (omega1, omega2, alpha1 - damping1, alpha2 - damping2)
  where
    delta = theta2 - theta1
    denom1 = (m1 + m2) * l1 - m2 * l1 * cos delta * cos delta
    denom2 = (l2 / l1) * denom1

    -- Angular accelerations
    alpha1 = (m2 * l1 * omega1 * omega1 * sin delta * cos delta
             + m2 * g * sin theta2 * cos delta
             + m2 * l2 * omega2 * omega2 * sin delta
             - (m1 + m2) * g * sin theta1) / denom1

    alpha2 = (- l2 * omega2 * omega2 * sin delta * cos delta
              + (m1 + m2) * g * sin theta1 * cos delta
              + (m1 + m2) * l1 * omega1 * omega1 * sin delta
              - (m1 + m2) * g * sin theta2) / denom2
    damping1 = frictionCoefficient * omega1 + dragCoefficient * omega1 * abs omega1
    damping2 = frictionCoefficient * omega2 + dragCoefficient * omega2 * abs omega2

-- Function to calculate the derivatives (right-hand side of the ODEs)
derivatives :: (Float, Float, Float, Float) -> (Float, Float, Float, Float)
derivatives (theta1, theta2, omega1, omega2) =
  (omega1, omega2, alpha1, alpha2)
  where
    delta = theta2 - theta1
    -- denom1 = (m1 + m2) * l1 - m2 * l1 * cos delta * cos delta
    -- denom2 = (l2 / l1) * denom1

    -- Angular accelerations
    -- Angular accelerations
    alpha1 = (-g * (2 * m1 + m2) * sin theta1
             - m2 * g * sin (theta1 - 2 * theta2)
             - 2 * sin delta * m2
               * (omega2^2 * l2 + omega1^2 * l1 * cos delta))
             / (l1 * (2 * m1 + m2 - m2 * cos (2 * delta)))
    

    alpha2 = (2 * sin delta
             * (omega1^2 * l1 * (m1 + m2)
             + g * (m1 + m2) * cos theta1
             + omega2^2 * l2 * m2 * cos delta))
             / (l2 * (2 * m1 + m2 - m2 * cos (2 * delta)))  

-- Air resistance and friction coefficients
dragCoefficient :: Float
dragCoefficient = 0.01  -- Example value

frictionCoefficient :: Float
frictionCoefficient = 0.01  -- Example value



-- Runge-Kutta 4th order step for double pendulum
rk4Step :: Float -> (Float, Float, Float, Float) -> (Float, Float, Float, Float)
rk4Step dt (theta1, theta2, omega1, omega2) =
  (theta1' , theta2', omega1', omega2')
  where
    -- Calculate k1
    (k1_theta1, k1_theta2, k1_omega1, k1_omega2) = derivatives' (theta1, theta2, omega1, omega2)

    -- Calculate k2
    (k2_theta1, k2_theta2, k2_omega1, k2_omega2) =
      derivatives' (theta1 + 0.5 * dt * k1_theta1, theta2 + 0.5 * dt * k1_theta2, 
                   omega1 + 0.5 * dt * k1_omega1, omega2 + 0.5 * dt * k1_omega2)

    -- Calculate k3
    (k3_theta1, k3_theta2, k3_omega1, k3_omega2) =
      derivatives' (theta1 + 0.5 * dt * k2_theta1, theta2 + 0.5 * dt * k2_theta2, 
                   omega1 + 0.5 * dt * k2_omega1, omega2 + 0.5 * dt * k2_omega2)

    -- Calculate k4
    (k4_theta1, k4_theta2, k4_omega1, k4_omega2) =
      derivatives' (theta1 + dt * k3_theta1, theta2 + dt * k3_theta2, 
                   omega1 + dt * k3_omega1, omega2 + dt * k3_omega2)

    -- Update the states using the RK4 formula
    theta1' = theta1 + (dt / 6) * (k1_theta1 + 2 * k2_theta1 + 2 * k3_theta1 + k4_theta1)
    theta2' = theta2 + (dt / 6) * (k1_theta2 + 2 * k2_theta2 + 2 * k3_theta2 + k4_theta2)
    omega1' = omega1 + (dt / 6) * (k1_omega1 + 2 * k2_omega1 + 2 * k3_omega1 + k4_omega1)
    omega2' = omega2 + (dt / 6) * (k1_omega2 + 2 * k2_omega2 + 2 * k3_omega2 + k4_omega2)
