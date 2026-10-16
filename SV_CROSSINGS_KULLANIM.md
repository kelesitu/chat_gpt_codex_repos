# MIMO tekil değer kesişimleri

MATLAB ile her frekansta büyükten küçüğe sıralanan tekil değerlerin seçilen dB seviyesindeki ilk, son ve tüm tespit edilen kesişim frekanslarını bulur.

## Kullanım

1. svCrossings.m ve runSvCrossings.m dosyalarını aynı klasöre indirin.
2. MATLAB Current Folder bölümünde bu klasörü açın.
3. runSvCrossings dosyasını çalıştırın. Örnek için Control System Toolbox gerekir.
4. Kendi sisteminiz için örnekteki L tanımını değiştirin.

Mevcut kompleks frekans cevaplarınızla doğrudan kullanım:

```matlab
% w: pozitif, artan frekans vektörü [rad/s]
% Hgrid: ny x nu x numel(w) kompleks frekans cevabı
% Hfun: tek frekansta kompleks cevap matrisi döndüren fonksiyon
[wFirst, wLast, wAll] = svCrossings(w, Hgrid, Hfun, -3);
```

Son parametre hedef seviyedir [dB]; verilmezse 0 dB kullanılır.
Sonuçlar rad/s birimindedir; Hz için 2*pi'ye bölün.
Kesişim yoksa ilk ve son sonuçları NaN olur.
wAll{k}, k'ncı en büyük tekil değerin tespit edilen tüm kesişimlerini içerir.
İsteğe bağlı dördüncü çıktı, −3 dB seviyesindeki son kesişimleri verir:

```matlab
[wFirst0, wLast0, wAll0, wLastMinus3] = svCrossings(w, Hgrid, Hfun, 0);
```

Sürekli zamanlı L için Hfun = @(wi) evalfr(L, 1i*wi).
Ayrık zamanlı L için Hfun = @(wi) freqresp(L, wi) kullanın;
w seçimini örnekleme zamanına ve Nyquist frekansına göre yapın.

## Yöntem ve sınırlar

Griddeki tekil değerler bir kez hesaplanır. İşaret değiştiren aralıklar,
logaritmik frekansta fzero ile hassaslaştırılır. TolX=1e-10 çözüm
toleransıdır; gerçek doğruluk modelin ve frekans cevabının sayısal
doğruluğuyla da sınırlıdır. Hgrid ve Hfun aynı sistemi temsil etmelidir.

İlk grid kesişiminden önce negatif eğimden pozitif eğime dönen çukur
adayları kaydedilir. En yakın adaydan başlayıp daha düşük frekanslara
doğru ilerlenir. Her çukur fminbnd ile hassaslaştırılır; çukurun iki
yanında işaret değişimi varsa fzero ile yeni kesişimler bulunur ve ilk
kesişim güncellenir. İlk bulunan yeni kesişimde arama durmaz.

Son grid kesişiminden sonra pozitif eğimden negatif eğime dönen tepe
adayları, en yakından başlayarak daha yüksek frekanslara doğru aynı
şekilde incelenir ve son kesişim güncellenir. Bu tarama −3 dB çıktısına
da uygulanır. Gridde hiç kesişim bulunmazsa bütün çukur ve tepe adayları
incelenir.

İlk ve son grid kesişimi arasındaki gizli çiftler için ayrıca ekstremum
taraması yapılmaz; wAll bütün matematiksel kökleri garanti etmez.
Gridde eğim değişimi izi bırakmayan dar salınımlar, aynı aday aralığındaki
birden fazla ekstremum ve örneklenmeyen teğet noktalar kaçabilir.
Sonlu grid ve yerel optimizasyon bütün kesişimleri garanti etmez.
Hedef seviyede düz bir bölüm varsa örneklenmiş noktalar kaydedilir;
bu bölüm izole bir kök değildir.

Tekil değerlerin sırası her frekansta büyükten küçüğedir; tekil vektör
dallarının devamlılığı takip edilmez.

runSvCrossings.m, analitik çözümü bilinen bir örnek ve doğrulama
assert ifadeleri içerir. Bu kayıt sırasında MATLAB ile çalıştırılmamıştır.
Çalıştırıldığında sonuçlar yerel svSonuclari.mat dosyasına kaydedilir.

## Yeni aday taraması kontrolleri

MATLAB'da `test_svCrossings` çalıştırılabilir; Control System Toolbox
gerektirmez. Kontroller, ilk/son kesişimin dışındaki gizli çiftleri,
hiç grid kesişimi olmayan çukur ve tepeleri, MIMO dallarını ve −3 dB
çıktısını kapsar. Bu güncelleme sırasında MATLAB bulunmadığından bu
dosya çalıştırılmadı; sayısal arama mantığının Python/SciPy karşılığı
analitik köklerle karşılaştırıldı.
