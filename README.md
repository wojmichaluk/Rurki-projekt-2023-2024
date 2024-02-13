# Rurki-projekt-2023-2024
### Repository for mini-project on RRiR (DE - Differential Equations) course at AGH University

Repozytorium zawiera moje rozwiązanie z użyciem MES - Metody Elementów Skończonych (**ang.** *FEM, Finite Element Method*) problemu 4.1 z pliku ```zadanie_obliczeniowe-KP-2.pdf```, którego "poprawiona" wersja jest przedstawiona w pliku ```corrected_task.png``` (w samym poleceniu nie ma stricte błędu, po prostu Podsiadło chciał, żeby problem do rozwiązania był nieco inny, a zapomniał zmienić plik). Rozwiązanie znajduje się w katalogu ```projekt``` - paczkę zip w dokładnie takiej postaci zamieściłem na UPEL i prezentowałem na zajęciach Podsiadle.
Projekt został oceniony na 50/50 punktów. Powinien być w pełni poprawny, choć zdaję sobie sprawę, że pewnie niektóre rzeczy dało się zrobić ładniej (np. te obliczenia na kartce w nienajlepszej jakości - nie jestem z nich dumny), no ale cóż :D. Przy prezentowaniu projektu prowadzący zadawał pytania m. in. o:
- wyprowadzenie sformułowania wariacyjnego (jak doszedłem do końcowej postaci macierzowej, gdzie skorzystałem z warunków brzegowych)
- sposób liczenia całki - w tym przypadku kwadratura Simpsona, a nie proponowana kwadratura Gaussa-Legandre'a
- obszar, na którym liczę całkę - w tym przypadku dwa "podprzedziały", na których funkcja podcałkowa jest niezerowa
- największy *n*, dla którego uzyskałem wyniki - w tym przypadku ~30_000, problemem jest pamięć a nie czas (przy nim odpalałem dla n = 10_000, kilkanaście sekund program liczył)
- jak szacuję złożoność programu - powiedziałem że ```O(n^2)```, zaakceptował - wydaje mi się, że tak jest
- poza tym sprawdzał jak wygladają wykresy, czy dla różnych *n* wyglądają mniej więcej tak samo

Mam nadzieję, że przedstawiona zawartość jest przydatna. Korzystajcie rozważnie!

UWAGA! Miałem temat 1 a nie 3, bo prowadzący przydzielał tematy po liście, nie wg algorytmu.
