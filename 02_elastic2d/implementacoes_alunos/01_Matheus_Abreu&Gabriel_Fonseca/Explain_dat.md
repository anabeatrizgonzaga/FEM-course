
Bloco 1: Parâmetros Globais de Controle:
  (1) número_de_nós (2) número_de_elementos (3) número_de_materiais (4) nós_por_elemento (5) gdl_por_nó (6) dimensão_do_problema
  
Bloco 2: Propriedades Restritas dos Materiais:
  (1) id_material (2) tipo_de_estado [1: EPD T3, 2: EPT T3, 3: EPD Q4, 4: EPT Q4]
  
Bloco 3: Propriedades Globais dos Materiais:
  (1) módulo_de_elasticidade (2) coeficiente_de_poisson (3) espessura [apenas de tipo_estado = 2]
   
Bloco 4: Coordenadas Nodais
  (1) id_nó (2) x_nó (3) y_nó (4) z_nó  
  
Bloco 5: Conectividade dos Elementos  
  (1) id_elemento (2) nó_1 (3) nó_2 (4) nó_3 (5) id_material
  
Bloco 6: Condições de Contorno (Restrições)  
  (1) id_nó  (2) gdl_1 (3) gdl_2 [0=livre; 1=fixo]
 
Bloco 7: Cargas Nodais 
  (1) id_nó  (2) valor_carga_x (3) valor_carga_y