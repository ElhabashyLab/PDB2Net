from py2cytoscape import cyrest
cy = cyrest.cyclient(port=1234)

# Teste, ob Cytoscape überhaupt ein leeres Netzwerk erstellen kann
print(cy.network.create())
