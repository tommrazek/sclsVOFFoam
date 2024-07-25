// solve Level-Set function as the re-initialization equation
   Info<< "solve the reinitialization equation"     
       << nl << endl;
   mixture.correctPsi0();
   psi == mixture.psi0();

   for (int corr=0; corr<int(mixture.epsilon().value()/mixture.deltaTau().value()); corr++)
   {
      psi = psi + mixture.psi0()/mag(mixture.psi0())*(double(1)-mag(fvc::grad(psi))*dimMeter)*mixture.deltaTau();
      psi.correctBoundaryConditions();
   }
