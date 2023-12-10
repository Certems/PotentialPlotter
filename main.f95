MODULE dataTypes
    !Types defined
    TYPE vortex
        REAL :: circulation=0.0
    END TYPE vortex
    TYPE source
        REAL :: strength=0.0
    END TYPE source
END MODULE dataTypes

program complexFluidPlotter
    !The program will;
    !1. Take a list of complex functions included in the plot
    !   1.1. Held in a list of compFunc types
    !2. Will place N particles within the domain (in grid formation)
    !   2.1. Held in a list of particle types
    !3. Calculate the motion of these particles over a time period
    !4. Record the motion in a list as it occurs (=> length proportional to the fixed timescale specified)
    !5. Output these position values as a text file (to be plotted in processing)    
    USE dataTypes
    implicit none

    !Useful variables initialised
    INTEGER, PARAMETER :: particleNumber = 5**2   !For simplicity, just deal with complete square arrangments => particle numbers square
    INTEGER, PARAMETER :: iterations     = 10
    TYPE(vortex), DIMENSION(1) :: potentialSet_vortex
    TYPE(source), DIMENSION(0) :: potentialSet_source
    COMPLEX, DIMENSION(particleNumber) :: particlePositions
    COMPLEX, DIMENSION(particleNumber, iterations) :: recordedPositions   !MAY HAVE WEIRD INDEXING, TRANSPOSE TO FIX
    INTEGER :: i

    !Computation begins
    PRINT *, "Computation START"
    !(1)
    PRINT *, "Step 1"
    CALL generate_potentialSet(potentialSet_vortex, potentialSet_source)
    !(2)
    PRINT *, "Step 2"
    CALL generate_particlePositions(particleNumber, particlePositions)
    !(3)
    PRINT *, "Step 3"
    DO i = 1, iterations
        PRINT *, "Step 4, Iteration ",i
        CALL record_particlePositions(i, particleNumber, iterations, particlePositions, recordedPositions) !(4)
        CALL calculate_iteration(particleNumber, particlePositions, potentialSet_vortex, potentialSet_source)
    END DO
    !(5)
    PRINT *, "Step 5"
    CALL output_data(recordedPositions)
    PRINT *, "Computation END"

CONTAINS

    SUBROUTINE generate_potentialSet(vortexSet, sourceSet)
        !Adds the list of potentials you wish to include to the list used to perform calculations
        !Important note; Ensure you change the size of the "potentialSet" array to MATCH the number of potentials you include here
        TYPE(vortex), DIMENSION(1), INTENT(INOUT) :: vortexSet
        TYPE(source), DIMENSION(1), INTENT(INOUT) :: sourceSet

        !Vortices
        TYPE(vortex) :: vortex_1
        vortex_1%circulation = 1.0
        vortexSet(1) = vortex_1

        !Sources
        !pass
    END SUBROUTINE generate_potentialSet

    SUBROUTINE generate_particlePositions(nParticles, positions)
        ! Generates the initial positions for the SQUARE arrangement of particles
        INTEGER, INTENT(IN) :: nParticles
        COMPLEX, DIMENSION(nParticles), INTENT(INOUT) :: positions
        REAL, PARAMETER :: scale = 1.0
        INTEGER :: j, i, index_j, index_i
        INTEGER :: squareLength
        COMPLEX :: initPos     !This still gives a deep copy as is strictly set = to in loop => considers it a new value, not overwritten perhaps
    
        ! Compute squareLength based on nParticles
        squareLength = CEILING(SQRT(REAL(nParticles)))
        initPos = (1.0, 4.0)
    
        ! Initialize positions
        index_j = 0
        DO j = -CEILING(squareLength/2.0), CEILING(squareLength/2.0)
            index_i = 0
            DO i = -CEILING(squareLength/2.0), CEILING(squareLength/2.0)
                initPos = COMPLEX(i*scale, j*scale)
                positions(index_j*squareLength +index_i) = initPos
                index_i = index_i +1
            END DO
            index_j = index_j +1
        END DO
    END SUBROUTINE generate_particlePositions

    SUBROUTINE calculate_iteration(nParticles, positions, vortexSet, sourceSet)
        !Calculates the motion of the z values in one cycle
        TYPE(vortex), DIMENSION(1), INTENT(IN) :: vortexSet
        TYPE(source), DIMENSION(0), INTENT(IN) :: sourceSet
        INTEGER, INTENT(IN) :: nParticles
        COMPLEX, DIMENSION(nParticles), INTENT(INOUT) :: positions
        INTEGER :: i,p,q
        COMPLEX :: newPosition
        DO i = 1,SIZE(positions)            !For each particle
            newPosition = positions(i)
            DO p = 1,SIZE(vortexSet)        !Sum contribution from all vortices
                newPosition = applyFunction_complex(newPosition, 0, vortexSet(p)%circulation)
            END DO
            DO q = 1,SIZE(sourceSet)        !Sum contribution from all sources
                newPosition = applyFunction_complex(newPosition, 1, sourceSet(q)%strength)
            END DO
            positions(i) = newPosition
        END DO
    END SUBROUTINE calculate_iteration

    SUBROUTINE record_particlePositions(currentIteration, nParticles, nIterations, positions, recorded_pos)
        !Takes current particle positions and stores them in an array (recordedPositions)
        INTEGER, INTENT(IN) :: currentIteration, nParticles, nIterations
        COMPLEX, DIMENSION(nParticles), INTENT(IN) :: positions
        COMPLEX, DIMENSION(nParticles, nIterations), INTENT(INOUT) :: recorded_pos
        INTEGER :: particleIndex
        DO particleIndex = 1,SIZE(positions)
            recorded_pos(particleIndex, currentIteration) = positions(particleIndex)
        END DO
    END SUBROUTINE record_particlePositions

    COMPLEX FUNCTION applyFunction_complex(z_0, functionPreset, functionMagnitude) RESULT(complexValue)
    !Takes a complex potential, applies it to the given point z_0, returns the new z
    COMPLEX :: z_0
    INTEGER :: functionPreset
    REAL :: functionMagnitude
    if(functionPreset == 0) then
        !Apply vortex
        complexValue = z_0 +EXP(-10.0/functionMagnitude)    !## DIFFERENT PSI FUNCTIONS WOULD RESULT IN DIFFERENT STEAMLINES, THIS ONLY EXPLORES A SINGLE STREAMLINE => LIKELY UNINTERESTING ##
    else if(functionPreset == 1) then
        !Apply line source
        complexValue = z_0 +EXP(10.0/functionMagnitude)     !## DIFFERENT PSI FUNCTIONS WOULD RESULT IN DIFFERENT STEAMLINES, THIS ONLY EXPLORES A SINGLE STREAMLINE => LIKELY UNINTERESTING ##
    END IF
    !if none of the above, then the function is invalid => ignore effect
    END FUNCTION applyFunction_complex

    SUBROUTINE output_data(recorded_pos)
        !Writes the recorded data to a textfile, to them be drawn in a separate program
        !   Writes position coordinates of each particle in rows with parts (RE, IM)
        !   Places a "-" on a line to signify the end of an iteration for all particles
        COMPLEX, DIMENSION(particleNumber, iterations) :: recorded_pos  !can read "particleNumber" & "iterations" from program above => not needed as arguement (same could be said the 2D array itself, but I am re-defining here to see if this is also valid or if an error is thrown)
        INTEGER :: j,i
        PRINT *, "Outputting Data..."
        OPEN(UNIT=1, FILE="outputData.txt", STATUS="replace", ACTION="write")
        DO j = 1,SIZE(recorded_pos, 2)
            DO i = 1,SIZE(recorded_pos, 1)
                WRITE(1, *) ",",REAL(recorded_pos(i,j)),",",AIMAG(recorded_pos(i,j))
            END DO
            WRITE(1, *) "-"
        END DO
        CLOSE(UNIT=1)
    END SUBROUTINE output_data

end program complexFluidPlotter