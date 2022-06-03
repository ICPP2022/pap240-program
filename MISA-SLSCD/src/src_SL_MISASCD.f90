!**********************************************************************************************************************
!>Main program (SLKMC algorithm)
!Function: it's used to simulate the process of radiation damage accumulation and subsequent annealing in metals
!          with the Synchronous Sublattice Kinetic Monte Carlo (SPKMC) algorithm.
!Date: 2019/12
!Author: CDD (chendandan@xs.ustb.edu.cn)
!Amendant Record: original code
!**********************************************************************************************************************
program SL_MISASCD
	use mod_constants
	use mod_structuretype
	use mod_globalvariables
	use mod_randdp
	implicit none
	include 'mpif.h'

	!<Used for communication
	type(defectCommHead) :: defCommList(3)
	logical, allocatable :: casCellList(:)          !<Used to mark the cell number that has implanted cascade
	integer :: arrayTemp(2,3)
	!Used for defects/reactions update
	integer, allocatable :: defUpdateStep(:,:)       !<defUpdateStep(SPECIES+6,numUpdates)
	integer :: numUpdates
	type(defect), pointer :: defectCurrent
	type(defectUpdate), pointer :: defUpdateList, defUpdate
	!<Record the choosed reaction
	type(reaction), pointer :: reacCH
	integer :: cellCH, lIndex
	type(defect), pointer :: defCH
	!<Record the runing time of the program
	double precision :: time1, time2, time3
	double precision :: commTime1, commTime2, runTime1, runTime2, postTime1, postTime2, postTime
	!<Used to open/close files
	character(len=12) :: form
	character(len=14) :: filename1, filename2, filename3, filename4
	!<Used for
	double precision :: numSendBuff(3), totalRecvBuff(3)
	!<Used for MPI
	integer :: status(MPI_STATUS_SIZE)
	!<Used to count simulations
	integer :: sim, outCount, vtkCount, nullSteps, totalCascades, cascadeCell, totCycle, sectorSteps(8)
	integer :: dim, sector, subCell, dir
	double precision :: totRateSector
	double precision :: tau_kmc, dt_kmc, t_kmc
	logical :: impCascadeToggle, irradToogle, alive
	!<Functions
	double precision, external :: compute_thresholdTime, SL_timeStepGenerate

	!test
	type(defect), pointer :: defTest
	type(reaction), pointer :: reacTest
	integer :: level, cellTest, jTest, countDef
	double precision :: countTest
	integer, external :: getLocalCell, find_level
	!*******************************************************

	interface
		subroutine SL_chooseReaction(sector, totRateSector, defCH, reacCH, cellCH)
			use mod_constants
			use mod_structuretype
			use mod_globalvariables
			use mod_randdp
			implicit none
			integer, intent(in) :: sector
			double precision, intent(in) :: totRateSector
			type(defect), pointer, intent(inout) :: defCH
			type(reaction), pointer, intent(inout) :: reacCH
			integer, intent(inout) :: cellCH
		end subroutine

		subroutine SL_updateDefectList(defCH,reacCH,cellCH,numUpdates,defUpdateStep,casCellList,defCommList)
			use mod_constants
			use mod_structuretype
			use mod_globalvariables
			use mod_randdp
			implicit none
			type(defect), pointer, intent(inout) :: defCH
			type(reaction), pointer, intent(in) :: reacCH
			integer, intent(in) :: cellCH, numUpdates
			integer, intent(inout), dimension(SPECIES+6,numUpdates) :: defUpdateStep
			logical, intent(inout), dimension(numMeshes) :: casCellList
			type(defectCommHead), intent(inout), dimension(3) :: defCommList
		end subroutine

		subroutine SL_updateReactionList(numUpdates,defUpdateStep)
			use mod_constants
			use mod_structuretype
			use mod_globalvariables
			implicit none
			integer, intent(in) :: numUpdates
			integer, intent(inout), dimension(SPECIES+6,numUpdates) :: defUpdateStep
		end subroutine

		subroutine SL_synDefectList(sector, defCommList, defUpdate, casCellList, arrayTemp)
			use mod_constants
			use mod_structuretype
			use mod_globalvariables
			implicit none
			integer, intent(in) :: sector
			type(defectCommHead), intent(inout), dimension(3) :: defCommList
			type(defectUpdate), pointer, intent(inout) :: defUpdate
			logical, intent(in), dimension(numMeshes) :: casCellList
			integer, intent(inout), dimension(2,3) :: arrayTemp
		end subroutine

		subroutine SL_synReactionList(defUpdateList)
			use mod_constants
			use mod_structuretype
			use mod_globalvariables
			implicit none
			type(defectUpdate), pointer, intent(in) :: defUpdateList
		end subroutine

		subroutine SL_synCascade(sector, casCellList, arrayTemp)
			use mod_constants
			use mod_structuretype
			use mod_globalvariables
			implicit none
			integer, intent(in) :: sector
			logical, intent(in), dimension(numMeshes) :: casCellList
			integer, intent(in), dimension(2,3) :: arrayTemp
		end subroutine
	end interface

	!<Start timing
	call cpu_time(time1)
	!*******************************************************
	!<Initialize MPI interface
	!*******************************************************
	ierr=0
	numtasks=0
	dims=0
	call MPI_INIT(ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
	call MPI_DIMS_CREATE(numtasks,3,dims, ierr)
	call MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,periods,.true.,comm,ierr)
	call MPI_CART_GET(comm,3,dims,periods,myProc%cartCoords,ierr)
	call MPI_CART_RANK(comm,myProc%cartCoords,myProc%taskid,ierr)
	call MPI_CART_SHIFT(comm,0,1,myProc%neighborProc(2),myProc%neighborProc(1),ierr) !x-direction
	call MPI_CART_SHIFT(comm,1,1,myProc%neighborProc(4),myProc%neighborProc(3),ierr) !y-direction
	call MPI_CART_SHIFT(comm,2,1,myProc%neighborProc(6),myProc%neighborProc(5),ierr) !z-direction

	!*******************************************************
	!<Read in parameters
	!*******************************************************
	call ReadInputs()				!Read in paramaters.txt, defect file, cascade file and mesh file.
	call SL_initializeMesh()
	allocate(casCellList(numMeshes))
	allocate(totalRateVol(numMeshes))		!<Create array of total rates in each mesh
	!<Initialize defectCommList. defectCommList is used to store defects that require communication
	do dim=1,3
		nullify(defCommList(dim)%commHead)
	end do

	if(myProc%taskid==0) then
		write(*,90) dims(1), dims(2), dims(3)
		90 format(' Proc division ', I3, I3, I3)
		write(*,100) myProc%taskid, numxLocal, numyLocal, numzLocal
		100 format(' Proc ', I3, ' numMeshes(x,y,z) ', I3, I3, I3)
		write(*,*)
	end if

	!******************************************************************************************************************!
	!******************************************************************************************************************!
	!*											MULTIPLE SIMULATION LOOP											  *!
	!*					Running multiple simulations. Initialize defects, reactions, boundary, etc.					  *!
	!******************************************************************************************************************!
	!******************************************************************************************************************!
	do sim=1, numSims
		!*******************************************************
		!<Initialize outputfiles
		!*******************************************************
		if(myProc%taskid==0) then
			select case(sim)
				case(1:9)
					write(form, '(I1)') sim
				case(10:99)
					write(form, '(I2)') sim
				case(100:999)
					write(form, '(I3)') sim
			end select
			if(totdatToggle=='yes') then	!<The same as MISA-SCD v1.0
				filename1(1:7)='totdat_'
				write(unit=filename1(8:10), fmt='(I3.3)') sim
				filename1(11:14)='.out'
				!write(filename1,*) "totdat_", trim(form),".out"
				open(TOTFILE, file=filename1, action='write', status='Unknown')
			end if
			if(defectToggle=='yes') then
				filename2(1:7)='defect_'
				write(unit=filename2(8:10), fmt='(I3.3)') sim
				filename2(11:14)='.out'
				open(DEFFILE, file=filename2, action='write', status='Unknown')
			end if
			if(stadatToggle=='yes') then
				filename3(1:7)='stadat_'
				write(unit=filename3(8:10), fmt='(I3.3)') sim
				filename3(11:14)='.out'
				open(STAFILE, file=filename3, action='write', status='Unknown')
			end if

			!<test
			!filename4 = 'oneCellDef.out'
			!open(15, file=filename4, action='write', status='Unknown')
		end if

		!*******************************************************
		!<Initialize random seed, defectList, reactionList, boundaryDefectList, totalRate.
		!*******************************************************
		call initializeRandomSeeds()	!<Initialize random seeds
		call compute_sinks()
		call initializeNumDefects()		!<Initialize point defects in the whole system
		call initializeDRL_Defect()		!<Build DRL for local area
		call initializeDRL_Ghost()		!<Build DRL for area
		call initializeDRL_Reaction()	!<Initialize reaction list and total rate of each mesh

		!*******************************************************
		!<Initialize post-processing parameters
		!*******************************************************
		irradToogle = .FALSE.
		if(totalDPA > 0d0 .AND. dpaRate > 0d0) then
			stopTime = totalDPA / dpaRate
			irradToogle = .TRUE.
			outputToggle = irradiationType
		else if(agingTime > 0d0) then
			stopTime = agingTime
			irradToogle = .TRUE.
			outputToggle = 'Aging'
		else
			stopTime = 0d0
			outputToggle = 'None'
		end if
		numImpEvents = 0
		numImpHe = 0
		numDamages = 0d0	!used to compute DPA
		time = 0d0			!<elapsed time
		step = 0
		nullSteps = 0
		totalCascades = 0
		outCount = 0
		vtkCount = 1

		!************************************************************************************************!
		!************************************************************************************************!
		!* 			                      		Damage Accumulation										*!
		!* During this process, Frenkerl Pairs or cascades are continuously implanted into system. 		*!
		!************************************************************************************************!
		!************************************************************************************************!
		totCycle=0
		allreduceTime = 0d0
		numCasTime = 0d0
		casCommTime = 0d0
		commTime11 = 0d0
		commTime12 = 0d0
		commTime21 = 0d0
		commTime22 = 0d0
		commCount_1 = 0
		commCount_2 = 0
		countEvent = 0

		countTest = 0d0
		runTime1 = MPI_WTIME()
		tau_kmc = INITAU
		do while(time < stopTime)
			!**********************************
			!<Sector LOOP
			!**********************************
			!tau_kmc = compute_thresholdTime()

			sectorSteps = 0
			do sector=1, 8
				casCellList=.FALSE.
				arrayTemp = 0

				!test test test
				totCycle = totCycle + 1

				!**********************************
				!<KMC LOOP
				!**********************************
				t_kmc = 0d0
				KMCLoop: do while(t_kmc < tau_kmc)
					cascadeCell = 0
					impCascadeToggle = .FALSE.

					dt_kmc = SL_timeStepGenerate(sector, totRateSector)	!totRateSector
					t_kmc = t_kmc + dt_kmc
					if(t_kmc > tau_kmc) then
						exit KMCLoop
					end if
					step = step + 1
					sectorSteps(sector) = sectorSteps(sector) + 1

					!<Choose a reaction in one mesh
					if(implantScheme=='explicit') then
						!if(time>=numImpEvents*(numDisplacedAtoms*atomVol)/(localVolume*dpaRate)) then
							!adding...
						!else
						!	call SL_chooseReaction(sector, totRateSector, defCH, reacCH, cellCH)
						!end if
					else if(implantScheme=='MonteCarlo') then
						call SL_chooseReaction(sector, totRateSector, defCH, reacCH, cellCH)
					end if

					if(.NOT. associated(reacCH)) then
						nullSteps = nullSteps + 1
						numUpdates = 0
					else
						if(reacCH%numReactants == -10) then	!0nd reaction: cascade
							numUpdates = 0
						else	!2nd reaction
							numUpdates = reacCH%numReactants + reacCH%numProducts	!= 2+1 or 2+2
						end if
						!<1:SPECIES: defectType; SPECIES+6: +/-1, cell, taskid, dir, neighbor, cascadeID
						allocate(defUpdateStep(SPECIES+6,numUpdates))

						!<Update defects according to reactions chosen.
						call SL_updateDefectList(defCH,reacCH,cellCH,numUpdates,defUpdateStep,&
								casCellList,defCommList)

						if(reacCH%numReactants == -10) then
							call resetReactionList(cellCH)
							cascadeCell = cellCH
							impCascadeToggle=.TRUE.
						else
							!<NOTE. This operation will change data in the memery that reacCH pointed to.
							call SL_updateReactionList(numUpdates,defUpdateStep)
						end if
						!<delete defUpdateStep(SPECIES+6, numUpdateDefects)
						deallocate(defUpdateStep)
					end if
				end do KMCLoop

				!**********************************
				!<Initialize defUpdateList. defUpdateList is used to store updated defects
				!**********************************
				allocate(defUpdateList)
				defUpdateList%cell = 0
				defUpdateList%proc = 0
				defUpdateList%dir = 0
				defUpdateList%neighbor = 0
				defUpdateList%cascadeID = 0
				nullify(defUpdateList%next)
				nullify(defUpdate)
				defUpdate => defUpdateList

				!**********************************
				!<Synchronization defects、reactions、cascaded efects
				!**********************************
				call SL_synDefectList(sector, defCommList, defUpdate, casCellList, arrayTemp)	!<create defUpdateList
			!	call MPI_BARRIER(comm)
				call SL_synReactionList(defUpdateList)					!<Update reactions according to defUpdateList
				if(irradiationType=='Cascade') then
					call SL_synCascade(sector, casCellList, arrayTemp)				!<Synchronization cascade defects and related diffusion reactions
				end if
				deallocate(defUpdateList)

			end do	!end sectors

			time = time + tau_kmc
			tau_kmc = compute_thresholdTime()
			!**********************************
			!<Post-processing process of damage accumulation
			!**********************************
			if(time >= stopTime/2.0d7*(2.0d0)**(outCount)) then
				numSendBuff(1) = dble(numImpEvents)
				numSendBuff(2) = dble(numImpHe)
				numSendBuff(3) = numDamages
				call MPI_REDUCE(numSendBuff, totalRecvBuff, 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0,comm, ierr)
				totImpEvents = totalRecvBuff(1)
				totImpHe = totalRecvBuff(2)
				totDamages = totalRecvBuff(3)
				call outputTotalDefects()	!<write totdat.out, defect.out, stadat.out

				runTime2 = MPI_WTIME()
				if(myProc%taskid==0) then
					if(irradiationType=='FrenkelPair') then
						DPA = dble(totImpEvents)/(systemVolume/(numDisplacedAtoms*atomVol))
					else if(irradiationType=='Cascade') then
						DPA = totDamages/(systemVolume/atomVol)
						DPA_low = dble(totImpEvents)/(systemVolume/(numDisplacedAtoms*atomVol))
					end if
					write(*,*) 'Time', time, 'totCycles', totCycle, 'CurrTau',tau_kmc, 'dt_kmc', dt_kmc
					if(irradiationType=='FrenkelPair') then
						write(*,*) 'DPA', DPA, 'FrenkelPair',totImpEvents, 'HeImpEvents', totImpHe
					else if(irradiationType=='Cascade') then
						write(*,*) 'DPA',DPA,'DPA_low',DPA_low,'Cascade',totImpEvents, 'HeImpEvents', totImpHe
					end if
					write(*,*)  'step', step,'runTime',runTime2-runTime1,&
							'commTime1',commTime11+commTime12,'commTime2',commTime21+commTime22
					write(*,*)
				end if
				outCount=outCount+1
			end if

			if(vtkToggle == 'yes') then
				if(time >= stopTime/10.0*vtkCount) then
					call outputVTK(vtkCount)
					vtkCount = vtkCount + 1
				end if
			end if

		end do	!end while(irradToogle .eqv. .true.)
		runTime2 = MPI_WTIME()
		!**********************************
		!<Final step: Post-processing process of damage accumulation
		!**********************************
		numSendBuff(1) = dble(numImpEvents)
		numSendBuff(2) = dble(numImpHe)
		numSendBuff(3) = numDamages
		call MPI_REDUCE(numSendBuff, totalRecvBuff, 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0,comm, ierr)
		totImpEvents = totalRecvBuff(1)
		totImpHe = totalRecvBuff(2)
		totDamages = totalRecvBuff(3)
		if(myProc%taskid == 0) then
			if(irradiationType=='FrenkelPair') then
				DPA = dble(totImpEvents)/(systemVolume/(numDisplacedAtoms*atomVol))
			else if(irradiationType=='Cascade') then
				DPA = totDamages/(systemVolume/atomVol)
				DPA_low = dble(totImpEvents)/(systemVolume/(numDisplacedAtoms*atomVol))
			end if
		end if
		call outputTotalDefects()	!<write totdat.out, defect.out, stadat.out

		!runTime2 = MPI_WTIME()
		if(myProc%taskid==0) then
			write(*,*) 'Final  step of damage accumulation'
			write(*,*) 'Time', time, 'totCycles', totCycle, 'CurrTau',tau_kmc, 'dt_kmc', dt_kmc
			if(irradiationType=='FrenkelPair') then
				write(*,*) 'DPA', DPA, 'FrenkelPair',totImpEvents, 'HeImpEvents', totImpHe
			else if(irradiationType=='Cascade') then
				write(*,*) 'DPA',DPA,'DPA_low',DPA_low,'Cascade',totImpEvents, 'HeImpEvents', totImpHe
			end if
			write(*,*)  'runTime',runTime2-runTime1,'commTime1',commTime11+commTime12,'commTime2',commTime21+commTime22
			write(*,*)
			write(*,*)
		end if


		!************************************************************************************************!
		!************************************************************************************************!
		!* 			                      		Annealing Loop											*!
		!* During this process, no Frenkerl Pairs or cascades are implanted into system. 				*!
		!************************************************************************************************!
		!************************************************************************************************!
		if(annealTime > 0d0) then
			outputToggle='Annealing'
			call initializeAnneal()
!			totalRate=totalRateCheck()

			time=0d0
			stopTime=annealTime
			annealStep=0
			outCount=0
			!<initialize anneal temperature
			annealIter=1
			temperature=annealTemperature
			casCellList=.FALSE.
			arrayTemp = 0

!			if(myProc%taskid==0) then
!				write(*,*)
!				write(*,*) 'Entering Annealing Phase'
!				if(totdatToggle=='yes') then
!					write(TOTFILE,*)
!					write(TOTFILE,*) 'Entering Annealing Phase'
!				end if
!				if(stadatToggle=='yes') then
!					write(STAFILE,*)
!					write(STAFILE,*) 'Entering Annealing Phase'
!				end if
!			end if

			do while(time < annealTime)
				!**********************************
				!<Sector LOOP
				!**********************************
				do sector=1, 8
					tau_kmc=compute_thresholdTime()

					!**********************************
					!<KMC LOOP
					!**********************************
					t_kmc = 0d0
					KMCLoop2: do while(t_kmc <= tau_kmc)
						if(time + t_kmc <= dble(annealIter)*annealTime/dble(annealSteps)) then

							!<Generate a timestep
							dt_kmc=SL_timeStepGenerate(sector, totRateSector)
							t_kmc= t_kmc + dt_kmc
							if(t_kmc > tau_kmc) then
								exit KMCLOOP2
							end if
							step=step+1
							annealStep=annealStep+1


							!<Choose a reaction in one mesh
							call SL_chooseReaction(sector, totRateSector, defCH, reacCH, cellCH)

							if(.NOT. associated(reacCH)) then
								nullSteps=nullSteps+1
								numUpdates=0
							else
								numUpdates=reacCH%numReactants+reacCH%numProducts
								allocate(defUpdateStep(SPECIES+6,numUpdates))
								!<Update defects according to reactions chosen.
								call SL_updateDefectList(defCH,reacCH,cellCH,numUpdates,&
										defUpdateStep,casCellList,defCommList)
								!<Update the corresponding reactions against defectCommList
								call SL_updateReactionList(numUpdates,defUpdateStep)
								!<delete defUpdateStep(SPECIES+6, numUpdateDefects)
								deallocate(defUpdateStep)
							end if
						else
							annealIter=annealIter+1
							if(annealIter <= annealSteps) then
								if(annealType=='mult') then
									temperature=annealTemperature*annealTempInc**dble(annealIter-1)
								else if(annealType=='add') then
									temperature=annealTemperature+annealTempInc*dble(annealIter-1)
								else if(annealType=='constant') then
									temperature=annealTemperature
								end if
							else
								exit KMCLoop2
							end if
						end if
					end do KMCLoop2	!<while( t_kmc <= tau_kmc)

					!**********************************
					!<Initialize defUpdateList. defUpdateList is used to store updated defects
					!**********************************
					allocate(defUpdateList)
					allocate(defUpdateList%defectType(SPECIES))
					defUpdateList%defectType=0
					nullify(defUpdateList%next)
					defUpdate=>defUpdateList
					!**********************************
					!<Synchronization defects and create defUpdateList
					!**********************************
					call SL_synDefectList(sector, defCommList, defUpdate, casCellList, arrayTemp)
					!**********************************
					!<Update reactions according to defUpdateList
					!**********************************
					call SL_synReactionList(defUpdateList)

					deallocate(defUpdateList%defectType)
					deallocate(defUpdateList)

				end do	!<sector=1, 8

				time = time + tau_kmc
				!**********************************
				!<Post-processing process of annealing
				!**********************************
				if(time >= annealTime/2.0d6*(2.0d0)**(outCount)) then
					call outputTotalDefects()	!<write totdat.out, defect.out, stadat.out
					!call outputTotalDefects_alternative()
					call cpu_time(time2)
					if(myProc%taskid==0) then
						write(*,*) 'Time', time, 'annealStep', annealStep
						write(*,*) outputToggle, 'runTime', time2-time1
						write(*,*)
					end if
					outCount=outCount+1
				end if

			end do	!<while(annealToggle .eqv. .TRUE.)

			!**********************************
			!<Final step: Post-processing process of annealing
			!**********************************
			call outputTotalDefects()	!<write totdat.out, defect.out, stadat.out
			call cpu_time(time2)
			if(myProc%taskid==0) then
				write(*,*) 'Final  step of annealing'
				write(*,*) 'Time', time, 'annealStep', annealStep
				write(*,*) outputToggle, 'runTime', time2-time1
			end if
		end if

		!***********************************************************************
		!Final step: deallocate defect lists and reaction lists
		!***********************************************************************
		call deallocateMyDRL()
		call deallocateMyGhostDRL()

		!<close files
		if(myProc%taskid==0) then
			write(*,*) 'Deallocating defects and reactions'
			if(totdatToggle=='yes') close(TOTFILE)
			if(defectToggle=='yes') close(DEFFILE)
			if(stadatToggle=='yes') close(STAFILE)
			!close(15)
		end if

	end do	!<end sim loop

	!***********************************************************************
	!<Clear memory
	!***********************************************************************
	deallocate(casCellList)
	deallocate(myMesh)
	do dir=1,6	!clear cells and loca in ghostMesh
		if(allocated(myGhost(dir)%cell)) then
			deallocate(myGhost(dir)%cell)
		end if
		if(allocated(myGhost(dir)%local)) then
			deallocate(myGhost(dir)%local)
		end if
	end do
	deallocate(totalRateVol)
	if(myProc%taskid==0) write(*,*) 'Deallocating memory: local meshes and ghost meshes'
	call deallocateInputs()
	call cpu_time(time2)
	if(myProc%taskid==0) then
		write(*,*) 'run time', runTime2-runTime1
		write(*,*) 'communication time', allreduceTime+numCasTime+casCommTime+commTime11+commTime12+commTime21+commTime22
		write(*,*) 'allresuce time', allreduceTime
		write(*,*) 'numCasTime',numCasTime,'casComm time', casCommTime
		write(*,*) 'pointComm time', commTime11+commTime12+commTime21+commTime22
		write(*,*)  'commTime11',commTime11,'commTime12',commTime12,'commTime21',commTime21,'commTime22',commTime22
		write(*,*) 'commCount_1', commCount_1
		write(*,*) 'number of boundary events', countEvent
	end if

	call MPI_FINALIZE(ierr)

end program
